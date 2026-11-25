#include "cubifier.hpp"
#include "constants.hpp"
#include "ids.hpp"
#include "shape_functions.hpp"
#include "signed_axis.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>

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

#include <stdexcept>
#include <vector>

#include <array>
#include <cassert>
#include <vector>

ElementMatrix Cubifier::computeHexaStiffness(const std::array<Vertex, 8> &X) {

  constexpr double a = 0.5773502691896257645;

  constexpr std::array<std::array<double, 3>, 8> signs = {{{{-1, -1, -1}},
                                                           {{1, -1, -1}},
                                                           {{1, 1, -1}},
                                                           {{-1, 1, -1}},
                                                           {{-1, -1, 1}},
                                                           {{1, -1, 1}},
                                                           {{1, 1, 1}},
                                                           {{-1, 1, 1}}}};

  ElementMatrix Ke = ElementMatrix::Zero();

  for (double xi : {-a, a}) {
    for (double eta : {-a, a}) {
      for (double zeta : {-a, a}) {

        // Shape function gradients in reference coordinates.
        Eigen::Matrix<double, 8, 3> dN;

        for (int i = 0; i < 8; ++i) {
          const auto &s = signs[i];

          const double sx = s[0];
          const double sy = s[1];
          const double sz = s[2];

          dN(i, 0) = sx * (1.0 + sy * eta) * (1.0 + sz * zeta) / 8.0;

          dN(i, 1) = sy * (1.0 + sx * xi) * (1.0 + sz * zeta) / 8.0;

          dN(i, 2) = sz * (1.0 + sx * xi) * (1.0 + sy * eta) / 8.0;
        }

        // Geometric Jacobian: dx/d(xi, eta, zeta).
        Eigen::Matrix3d J = Eigen::Matrix3d::Zero();

        for (int i = 0; i < 8; ++i) {
          J.noalias() += X[i] * dN.row(i);
        }

        const double detJ = J.determinant();

        if (!std::isfinite(detJ) || detJ <= 1e-12) {
          throw std::runtime_error("Invalid HEX8 Jacobian");
        }

        // Physical gradients:
        // gradN = dN * J^{-1}.
        //
        // Solve J^T * gradN^T = dN^T instead of
        // explicitly computing J.inverse().

        Eigen::Matrix<double, 8, 3> gradN =
            J.transpose().fullPivLu().solve(dN.transpose()).transpose();

        // Integrate grad(N_i) . grad(N_j).
        Ke.noalias() += gradN * gradN.transpose() * detJ;
      }
    }
  }

  return Ke;
}

void Cubifier::build() {
  vertices.clear();
  hexes.clear();
  constrainedFlag.clear();

  const auto start = region.startSlices;
  const auto diff = region.endSlices - start;

  using NodeIdsTensor = Eigen::Tensor<NodeId, 3>;

  NodeIdsTensor nodeIds(diff.x() + 1, diff.y() + 1, diff.z() + 1);

  nodeIds.setConstant(NO_NODE);

  // ------------------------------------------------------------
  // 1. Create vertices and HEX8 elements
  // ------------------------------------------------------------

  const std::array<Eigen::Vector3i, 8> offsets = {{{0, 0, 0},
                                                   {1, 0, 0},
                                                   {1, 1, 0},
                                                   {0, 1, 0},
                                                   {0, 0, 1},
                                                   {1, 0, 1},
                                                   {1, 1, 1},
                                                   {0, 1, 1}}};

  for (const auto &voxel : region.voxels) {

    Hex hex;

    for (int i = 0; i < 8; ++i) {

      const auto pos = voxel + offsets[i];

      NodeId &id = nodeIds(pos.x() - start.x(), pos.y() - start.y(),
                           pos.z() - start.z());

      if (id == NO_NODE) {

        id = static_cast<NodeId>(vertices.size());

        // Keep original physical coordinates.
        // Do not normalize them here.
        vertices.push_back(Vertex(static_cast<double>(pos.x()),
                                  static_cast<double>(pos.y()),
                                  static_cast<double>(pos.z())));
      }

      hex[i] = id;
    }

    hexes.push_back(hex);
  }

  const NodeId n = static_cast<NodeId>(vertices.size());

  if (n == 0) {
    throw std::runtime_error("No vertices were generated");
  }

  // ------------------------------------------------------------
  // 2. Collect boundary nodes
  // ------------------------------------------------------------

  std::array<std::array<std::vector<NodeId>, 2>, 3> boundaryNodes;

  auto getBoundaryNodes = [&](Axis a,
                              IntegerSign sign) -> std::vector<NodeId> & {
    const int axis = static_cast<int>(a);
    const int side = (sign == POSITIVE) ? 0 : 1;

    return boundaryNodes[axis][side];
  };

  // ------------------------------------------------------------
  // 2a. Boundary nodes identified by the surface mesh
  // ------------------------------------------------------------

  const auto &quadIds = surfaceMesh.regionToQuadIds[region.id];

  std::array<IntegerSign, 3> missingAxis = {UNDEFINED, UNDEFINED, UNDEFINED};

  for (const auto quadId : quadIds) {

    const Axis a = surfaceMesh.quadToAxisConstraint[quadId];

    const IntegerSign sign = region.facesNormals[a];

    // The opposite face may need to be added from the
    // region's bounding box.
    missingAxis[a] = (sign == POSITIVE) ? NEGATIVE : POSITIVE;

    for (const auto surfaceNodeId : surfaceMesh.quads[quadId]) {

      auto nodePos = surfaceMesh.nodes[surfaceNodeId] - start;

      const NodeId id = nodeIds(nodePos.x(), nodePos.y(), nodePos.z());

      if (id == NO_NODE) {
        throw std::runtime_error(
            "Surface node does not belong to the voxel region");
      }

      getBoundaryNodes(a, sign).push_back(id);
    }
  }

  // ------------------------------------------------------------
  // 2b. Add missing bounding-box faces
  // ------------------------------------------------------------

  for (Axis a = 0; a < DIMENSION; ++a) {

    for (IntegerSign sign : {POSITIVE, NEGATIVE}) {

      if (missingAxis[a] != sign && missingAxis[a] != UNDEFINED)
        continue;

      Voxel min = Eigen::Vector3i::Zero();
      Voxel max = diff;

      if (sign == POSITIVE) {
        min[a] = max[a] - 1;
      } else {
        max[a] = 1;
      }

      for (int x = min.x(); x <= max.x(); ++x) {
        for (int y = min.y(); y <= max.y(); ++y) {
          for (int z = min.z(); z <= max.z(); ++z) {

            const NodeId id = nodeIds(x, y, z);

            if (id != NO_NODE) {
              getBoundaryNodes(a, sign).push_back(id);
            }
          }
        }
      }
    }
  }

  // ------------------------------------------------------------
  // 3. Build coordinate-wise Dirichlet constraints
  // ------------------------------------------------------------

  Eigen::Matrix<bool, Eigen::Dynamic, 3> constrained =
      Eigen::Matrix<bool, Eigen::Dynamic, 3>::Constant(n, 3, false);

  Eigen::MatrixXd boundaryValues = Eigen::MatrixXd::Zero(n, 3);

  for (Axis a = 0; a < DIMENSION; ++a) {

    const int coordinate = static_cast<int>(a);

    for (IntegerSign sign : {POSITIVE, NEGATIVE}) {

      const auto &nodes = getBoundaryNodes(a, sign);

      const double target = (sign == POSITIVE) ? 1.0 : 0.0;

      for (NodeId id : nodes) {

        if (id >= n) {
          throw std::runtime_error("Invalid boundary node ID");
        }

        if (constrained(id, coordinate)) {

          const double error =
              std::abs(boundaryValues(id, coordinate) - target);

          if (error > 1e-8) {
            throw std::runtime_error("Inconsistent boundary constraints");
          }

        } else {

          constrained(id, coordinate) = true;
          boundaryValues(id, coordinate) = target;
        }
      }
    }
  }

  // Keep this only if constrainedFlag is used elsewhere.
  // It now indicates whether a node has at least one
  // constrained coordinate.
  constrainedFlag.resize(n);

  for (NodeId i = 0; i < n; ++i) {
    // constrainedFlag[i] =
    //     constrained(i, 0) || constrained(i, 1) || constrained(i, 2);
    constrainedFlag[i] = constrained(i, 2);
  }

  // ------------------------------------------------------------
  // 4. Assemble global stiffness matrix
  // ------------------------------------------------------------

  std::vector<Eigen::Triplet<double>> triplets;

  triplets.reserve(hexes.size() * 64);

  for (const auto &hex : hexes) {

    std::array<Vertex, 8> elementVertices;

    for (int i = 0; i < 8; ++i) {
      elementVertices[i] = vertices[hex[i]];
    }

    const ElementMatrix Ke = computeHexaStiffness(elementVertices);

    for (int i = 0; i < 8; ++i) {
      for (int j = 0; j < 8; ++j) {

        triplets.emplace_back(hex[i], hex[j], Ke(i, j));
      }
    }
  }

  Eigen::SparseMatrix<double> K(n, n);

  K.setFromTriplets(triplets.begin(), triplets.end());

  K.makeCompressed();

  // ------------------------------------------------------------
  // 5. Solve independently for u, v and w
  // ------------------------------------------------------------

  Eigen::MatrixXd mappedVertices = Eigen::MatrixXd::Zero(n, 3);

  for (int coordinate = 0; coordinate < 3; ++coordinate) {

    std::vector<NodeId> interior;
    std::vector<int> globalToInterior(n, -1);

    // Find vertices whose current coordinate is unknown.
    for (NodeId i = 0; i < n; ++i) {

      if (!constrained(i, coordinate)) {

        globalToInterior[i] = static_cast<int>(interior.size());

        interior.push_back(i);
      }
    }

    const int ni = static_cast<int>(interior.size());

    // Store all prescribed values for this coordinate.
    for (NodeId i = 0; i < n; ++i) {

      if (constrained(i, coordinate)) {
        mappedVertices(i, coordinate) = boundaryValues(i, coordinate);
      }
    }

    if (ni == 0) {
      continue;
    }

    // --------------------------------------------------------
    // Build reduced system:
    //
    // K_II * x_I = -K_IB * x_B
    // --------------------------------------------------------

    std::vector<Eigen::Triplet<double>> reducedTriplets;

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(ni);

    for (int i = 0; i < ni; ++i) {

      const NodeId globalI = interior[i];

      // Eigen's default SparseMatrix is column-major.
      // InnerIterator(K, globalI) iterates column globalI.
      // Therefore, the varying index is it.row().
      for (Eigen::SparseMatrix<double>::InnerIterator it(K, globalI); it;
           ++it) {

        const NodeId globalJ = it.row();

        const double value = it.value();

        if (constrained(globalJ, coordinate)) {

          rhs[i] -= value * boundaryValues(globalJ, coordinate);

        } else {

          reducedTriplets.emplace_back(i, globalToInterior[globalJ], value);
        }
      }
    }

    Eigen::SparseMatrix<double> A(ni, ni);

    A.setFromTriplets(reducedTriplets.begin(), reducedTriplets.end());

    A.makeCompressed();

    // --------------------------------------------------------
    // Solve
    // --------------------------------------------------------

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    solver.compute(A);

    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Failed to factorize harmonic system");
    }

    const Eigen::VectorXd solution = solver.solve(rhs);

    if (solver.info() != Eigen::Success || !solution.allFinite()) {
      throw std::runtime_error("Failed to solve harmonic system");
    }

    // --------------------------------------------------------
    // Store this coordinate
    // --------------------------------------------------------

    for (int i = 0; i < ni; ++i) {

      mappedVertices(interior[i], coordinate) = solution[i];
    }
  }

  // ------------------------------------------------------------
  // 6. Replace original coordinates with mapped coordinates
  // ------------------------------------------------------------

  for (NodeId i = 0; i < n; ++i) {

    vertices[i] = mappedVertices.row(i).transpose();
  }
}

void Cubifier::exportMesh(std::string filename) const {
  vtkNew<vtkPoints> points;
  points->SetDataTypeToDouble();
  points->SetNumberOfPoints(vertices.size());

  // Add vertices.
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    const auto &v = vertices[i];
    points->SetPoint(i, v.x(), v.y(), v.z());
  }

  // Create unstructured grid.
  vtkNew<vtkUnstructuredGrid> grid;
  grid->SetPoints(points);

  // Add hexahedral cells.
  vtkNew<vtkIdList> pointIds;

  for (const auto &hex : hexes) {
    pointIds->Reset();

    for (const auto nodeId : hex) {
      pointIds->InsertNextId(static_cast<vtkIdType>(nodeId));
    }

    grid->InsertNextCell(VTK_HEXAHEDRON, pointIds);
  }

  // Add constrainedFlag as point data.
  if (constrainedFlag.size() != vertices.size()) {
    throw std::runtime_error(
        "constrainedFlag size does not match vertices size");
  }

  vtkNew<vtkUnsignedCharArray> constrainedArray;
  constrainedArray->SetName("constrainedFlag");
  constrainedArray->SetNumberOfComponents(1);
  constrainedArray->SetNumberOfTuples(vertices.size());

  size_t nb = 0;
  for (std::size_t i = 0; i < constrainedFlag.size(); ++i) {
    if (constrainedFlag[i])
      nb++;
    constrainedArray->SetValue(static_cast<vtkIdType>(i),
                               constrainedFlag[i] ? 1 : 0);
  }

  grid->GetPointData()->AddArray(constrainedArray);

  // Write VTU file.
  vtkNew<vtkXMLUnstructuredGridWriter> writer;
  writer->SetFileName(filename.c_str());
  writer->SetInputData(grid);
  writer->SetDataModeToBinary();

  if (writer->Write() != 1) {
    throw std::runtime_error("Failed to write mesh to " + filename);
  }
}
