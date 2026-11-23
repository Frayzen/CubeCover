#include "cubifier.hpp"
#include "ids.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <stdexcept>
#include <string>

#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <algorithm>
#include <array>
#include <set>
#include <stdexcept>
#include <vector>

void Cubifier::build() {
  using SparseMatrix = Eigen::SparseMatrix<double>;
  using Triplet = Eigen::Triplet<double>;

  // ------------------------------------------------------------
  // 1. Region dimensions
  // ------------------------------------------------------------

  const Voxel start = region.startSlices;
  const Voxel diff = region.endSlices - start;

  const int nx = diff.x();
  const int ny = diff.y();
  const int nz = diff.z();

  if (nx <= 0 || ny <= 0 || nz <= 0)
    throw std::runtime_error("Invalid region dimensions.");

  const int nodeNx = nx + 1;
  const int nodeNy = ny + 1;
  const int nodeNz = nz + 1;

  const std::size_t totalGridNodes =
      static_cast<std::size_t>(nodeNx) * nodeNy * nodeNz;

  // ------------------------------------------------------------
  // 2. Voxel occupancy
  // ------------------------------------------------------------

  Eigen::Tensor<bool, DIMENSION> voxelFlag(nx, ny, nz);
  Eigen::Tensor<VoxelId, DIMENSION> voxelId(nx, ny, nz);

  voxelFlag.setConstant(false);

  VoxelId curId = 0;

  for (Voxel v : region) {
    const Voxel rel = v - start;

    voxelFlag(rel.x(), rel.y(), rel.z()) = true;
    voxelId(rel.x(), rel.y(), rel.z()) = curId++;
  }

  // ------------------------------------------------------------
  // 3. Grid nodes and graph edges
  // ------------------------------------------------------------

  vertices.clear();
  nodeGridCoordinates.clear();

  std::vector<int> nodeLookup(totalGridNodes, -1);

  std::vector<std::array<double, 3>> fixedValues;
  std::vector<std::array<bool, 3>> fixedMask;

  auto flattenNode = [&](const Eigen::Vector3i &p) -> std::size_t {
    return (static_cast<std::size_t>(p.x()) * nodeNy + p.y()) * nodeNz + p.z();
  };

  auto getNodeId = [&](const Eigen::Vector3i &p) -> int {
    const std::size_t index = flattenNode(p);

    if (nodeLookup[index] != -1)
      return nodeLookup[index];

    const int id = static_cast<int>(nodeGridCoordinates.size());

    nodeLookup[index] = id;

    nodeGridCoordinates.push_back(p);
    fixedValues.push_back({0.0, 0.0, 0.0});
    fixedMask.push_back({false, false, false});

    return id;
  };

  // Each edge is stored only once.
  std::set<std::pair<int, int>> edges;

  auto addEdge = [&](int a, int b) {
    if (a > b)
      std::swap(a, b);

    edges.emplace(a, b);
  };

  // Each occupied voxel contributes its 8 corners and 12 edges.
  for (int x = 0; x < nx; ++x) {
    for (int y = 0; y < ny; ++y) {
      for (int z = 0; z < nz; ++z) {
        if (!voxelFlag(x, y, z))
          continue;

        const Eigen::Vector3i base(x, y, z);

        int corners[8];

        for (int b = 0; b < 8; ++b) {
          Eigen::Vector3i p = base;

          p.x() += (b >> 0) & 1;
          p.y() += (b >> 1) & 1;
          p.z() += (b >> 2) & 1;

          corners[b] = getNodeId(p);
        }

        // Connect corners differing along exactly one axis.
        for (int b = 0; b < 8; ++b) {
          for (int a = 0; a < 3; ++a) {
            if ((b & (1 << a)) == 0) {
              addEdge(corners[b], corners[b | (1 << a)]);
            }
          }
        }
      }
    }
  }

  const int numNodes = static_cast<int>(nodeGridCoordinates.size());

  if (numNodes == 0)
    throw std::runtime_error("Region contains no voxels.");

  // Build adjacency.
  std::vector<std::vector<int>> adjacency(numNodes);

  for (const auto &[a, b] : edges) {
    adjacency[a].push_back(b);
    adjacency[b].push_back(a);
  }

  // ------------------------------------------------------------
  // 4. Boundary constraints
  // ------------------------------------------------------------

  // Pair indexing:
  // signIndex = 0 -> NEGATIVE
  // signIndex = 1 -> POSITIVE

  std::array<std::array<bool, 2>, 3> encountered{};

  auto getSignIndex = [](IntegerSign sign) -> int {
    if (sign == NEGATIVE)
      return 0;

    if (sign == POSITIVE)
      return 1;

    throw std::runtime_error(
        "Invalid surface normal sign: expected POSITIVE or NEGATIVE.");
  };

  auto setConstraint = [&](int nodeId, int axis, double value) {
    if (fixedMask[nodeId][axis]) {
      if (std::abs(fixedValues[nodeId][axis] - value) > 1e-10) {
        throw std::runtime_error(
            "Conflicting boundary constraints on grid node.");
      }

      return;
    }

    fixedMask[nodeId][axis] = true;
    fixedValues[nodeId][axis] = value;
  };

  // Constrain the four nodes belonging to a voxel face.
  auto constrainVoxelFace = [&](const Eigen::Vector3i &voxel, int axis,
                                int signIndex) {
    Eigen::Vector3i faceOrigin = voxel;

    // NEGATIVE: lower face of voxel.
    // POSITIVE: upper face of voxel.
    if (signIndex == 1)
      faceOrigin[axis] += 1;

    const double target = static_cast<double>(signIndex);

    int otherAxes[2];
    int n = 0;

    for (int a = 0; a < 3; ++a) {
      if (a != axis)
        otherAxes[n++] = a;
    }

    for (int u = 0; u <= 1; ++u) {
      for (int v = 0; v <= 1; ++v) {
        Eigen::Vector3i p = faceOrigin;

        p[otherAxes[0]] += u;
        p[otherAxes[1]] += v;

        const int id = getNodeId(p);

        setConstraint(id, axis, target);
      }
    }
  };

  // ------------------------------------------------------------
  // 4.1 Apply existing surface constraints
  // ------------------------------------------------------------

  for (QuadId qId : surfaceMesh.regionToQuads[region.id]) {
    const Voxel v = surfaceMesh.quadToVoxel[qId];

    const Axis axis = surfaceMesh.quadToAxisConstraint[qId];
    const IntegerSign sign = region.facesNormals[axis];

    if (axis < 0 || axis >= 3)
      throw std::runtime_error("Invalid surface normal axis.");

    const int signIndex = getSignIndex(sign);

    const Voxel rel = v - start;

    if (rel.x() < 0 || rel.x() >= nx || rel.y() < 0 || rel.y() >= ny ||
        rel.z() < 0 || rel.z() >= nz) {
      throw std::runtime_error(
          "Surface quad references a voxel outside the region.");
    }

    if (!voxelFlag(rel.x(), rel.y(), rel.z()))
      throw std::runtime_error("Surface quad references an unoccupied voxel.");

    encountered[axis][signIndex] = true;

    constrainVoxelFace(Eigen::Vector3i(rel.x(), rel.y(), rel.z()), axis,
                       signIndex);
  }

  // ------------------------------------------------------------
  // 4.2 Complete missing boundary constraints
  // ------------------------------------------------------------

  for (int axis = 0; axis < 3; ++axis) {
    for (int signIndex = 0; signIndex < 2; ++signIndex) {
      if (encountered[axis][signIndex])
        continue;

      const int boundaryIndex = (signIndex == 0) ? 0
                                                 : ((axis == 0)   ? nx - 1
                                                    : (axis == 1) ? ny - 1
                                                                  : nz - 1);

      for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
          for (int z = 0; z < nz; ++z) {
            if (!voxelFlag(x, y, z))
              continue;

            const Eigen::Vector3i voxel(x, y, z);

            if (voxel[axis] != boundaryIndex)
              continue;

            constrainVoxelFace(voxel, axis, signIndex);
          }
        }
      }
    }
  }

  // ------------------------------------------------------------
  // 5. Solve the three harmonic coordinate systems
  // ------------------------------------------------------------

  vertices.assign(numNodes, Eigen::Vector3d::Zero());

  for (int coordinate = 0; coordinate < 3; ++coordinate) {
    std::vector<int> interiorNodes;
    std::vector<int> interiorIndex(numNodes, -1);

    for (int i = 0; i < numNodes; ++i) {
      if (!fixedMask[i][coordinate]) {
        interiorIndex[i] = static_cast<int>(interiorNodes.size());

        interiorNodes.push_back(i);
      } else {
        vertices[i][coordinate] = fixedValues[i][coordinate];
      }
    }

    const int numInterior = static_cast<int>(interiorNodes.size());
    std::cout << "GOT " << numInterior << " FOR AXIS " << coordinate
              << std::endl;

    if (numInterior == 0)
      continue;

    std::vector<Triplet> triplets;
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(numInterior);

    // Assemble the reduced Laplacian:
    //
    // L_II * x_I = -L_IC * x_C

    for (int row = 0; row < numInterior; ++row) {
      const int node = interiorNodes[row];

      const auto &neighbors = adjacency[node];

      const double degree = static_cast<double>(neighbors.size());

      if (degree == 0.0)
        throw std::runtime_error("Isolated grid node encountered.");

      // Diagonal term.
      triplets.emplace_back(row, row, degree);

      for (int neighbor : neighbors) {
        if (fixedMask[neighbor][coordinate]) {
          // Move fixed values to the right-hand side.
          rhs[row] += fixedValues[neighbor][coordinate];
        } else {
          const int col = interiorIndex[neighbor];

          triplets.emplace_back(row, col, -1.0);
        }
      }
    }

    SparseMatrix A(numInterior, numInterior);

    A.setFromTriplets(triplets.begin(), triplets.end());

    A.makeCompressed();

    Eigen::SimplicialLDLT<SparseMatrix> solver;

    solver.compute(A);

    if (solver.info() != Eigen::Success)
      throw std::runtime_error("Failed to factorize the Laplacian matrix.");

    const Eigen::VectorXd solution = solver.solve(rhs);

    if (solver.info() != Eigen::Success || !solution.allFinite()) {
      throw std::runtime_error("Failed to solve the Laplacian system.");
    }

    for (int i = 0; i < numInterior; ++i) {
      vertices[interiorNodes[i]][coordinate] = solution[i];
    }
  }
}

void Cubifier::exportMesh(std::string filename) const {

  // --------------------------------------------------------
  // Validate node data
  // --------------------------------------------------------

  if (!nodeGridCoordinates.empty() &&
      nodeGridCoordinates.size() != vertices.size()) {
    throw std::runtime_error("Cubifier: nodeGridCoordinates size mismatch.");
  }

  // --------------------------------------------------------
  // Points
  // --------------------------------------------------------

  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(static_cast<vtkIdType>(vertices.size()));

  for (NodeId i = 0; i < vertices.size(); ++i) {

    const Vertex &v = vertices[i];

    points->SetPoint(static_cast<vtkIdType>(i), v.x(), v.y(), v.z());
  }

  // --------------------------------------------------------
  // Unstructured grid
  // --------------------------------------------------------

  vtkNew<vtkUnstructuredGrid> grid;
  grid->SetPoints(points);

  // --------------------------------------------------------
  // Hexahedra
  // --------------------------------------------------------

  vtkNew<vtkHexahedron> hex;

  for (const Hex &cell : hexes) {

    if (cell.size() != 8) {
      throw std::runtime_error("Cubifier: a hexahedron must have 8 vertices.");
    }

    for (int i = 0; i < 8; ++i) {

      const NodeId vertexId = cell[i];

      if (vertexId >= vertices.size()) {
        throw std::runtime_error("Cubifier: invalid hexahedron vertex index.");
      }

      hex->GetPointIds()->SetId(i, static_cast<vtkIdType>(vertexId));
    }

    grid->InsertNextCell(VTK_HEXAHEDRON, hex->GetPointIds());
  }

  // --------------------------------------------------------
  // Node positions (point data)
  // --------------------------------------------------------

  if (!vertices.empty()) {

    vtkNew<vtkDoubleArray> positions;
    positions->SetName("NodePosition");
    positions->SetNumberOfComponents(3);
    positions->SetNumberOfTuples(static_cast<vtkIdType>(vertices.size()));

    for (vtkIdType i = 0; i < vertices.size(); ++i) {

      const Eigen::Vector3d &p = vertices[i];

      positions->SetTuple3(i, p.x(), p.y(), p.z());
    }

    grid->GetPointData()->AddArray(positions);
  }

  // --------------------------------------------------------
  // Node grid coordinates (point data)
  // --------------------------------------------------------

  if (!nodeGridCoordinates.empty()) {

    vtkNew<vtkIntArray> coordinates;
    coordinates->SetName("NodeGridCoordinates");
    coordinates->SetNumberOfComponents(3);
    coordinates->SetNumberOfTuples(
        static_cast<vtkIdType>(nodeGridCoordinates.size()));

    for (vtkIdType i = 0; i < nodeGridCoordinates.size(); ++i) {

      const Eigen::Vector3i &p = nodeGridCoordinates[i];

      coordinates->SetTuple3(i, p.x(), p.y(), p.z());
    }

    grid->GetPointData()->AddArray(coordinates);
  }

  // --------------------------------------------------------
  // Write VTU
  // --------------------------------------------------------

  vtkNew<vtkXMLUnstructuredGridWriter> writer;

  writer->SetFileName(filename.c_str());
  writer->SetInputData(grid);
  writer->SetDataModeToBinary();

  if (writer->Write() != 1) {
    throw std::runtime_error("Cubifier: failed to write VTU file: " + filename);
  }
}
