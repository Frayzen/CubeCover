#include "cubifier.hpp"
#include "constants.hpp"
#include "ids.hpp"
#include "signed_axis.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <stdexcept>
#include <string>

#include <unordered_set>
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
#include <unordered_set>
#include <vector>

void Cubifier::applyConstraints(const Eigen::Tensor<NodeId, DIMENSION> &nodeIds,
                                const Voxel &start,
                                std::vector<bool> &constrainedFlag, bool last) {

  std::array<Axis, DIMENSION> internalFaceSign = {UNDEFINED, UNDEFINED,
                                                  UNDEFINED};

  // Apply constraints from the surface mesh.
  for (QuadId qId : surfaceMesh.regionToQuads[region.id]) {
    const Axis axis = surfaceMesh.quadToNormal[qId];
    const double value = getIntegerSignValue(region.facesNormals[axis]);

    const Quad &q = surfaceMesh.quads[qId];

    for (NodeId surfaceNodeId : q) {
      const Voxel v = surfaceMesh.nodes[surfaceNodeId];
      const Voxel relPos = v - start;

      const NodeId systemNodeId = nodeIds(relPos.x(), relPos.y(), relPos.z());

      assert(systemNodeId != NO_NODE);
      assert(systemNodeId < vertices.size());

      constrainedFlag[systemNodeId] = !last;
      vertices[systemNodeId][axis] = value;
    }

    internalFaceSign[axis] =
        region.facesNormals[axis] == POSITIVE ? NEGATIVE : POSITIVE;
  }

  for (Axis axis = 0; axis < DIMENSION; ++axis) {
    for (IntegerSign sign : {POSITIVE, NEGATIVE}) {
      if (internalFaceSign[axis] != UNDEFINED && internalFaceSign[axis] != sign)
        continue;
      const double value = getIntegerSignValue(sign);

      Voxel min = Voxel::Zero();
      Voxel max = Voxel(nodeIds.dimension(0), nodeIds.dimension(1),
                        nodeIds.dimension(2));

      if (sign == POSITIVE) {
        min[axis] = nodeIds.dimension(axis) - 1;
      } else {
        max[axis] = 1;
      }

      for (int x = min.x(); x < max.x(); ++x) {
        for (int y = min.y(); y < max.y(); ++y) {
          for (int z = min.z(); z < max.z(); ++z) {
            const NodeId id = nodeIds(x, y, z);
            if (id == NO_NODE)
              continue;
            vertices[id][axis] = value;
            constrainedFlag[id] = true;
          }
        }
      }
    }
  }
}

void Cubifier::build() {
  vertices.clear();
  hexes.clear();
  constrainedFlag.clear();

  const Voxel start = region.startSlices;
  const Voxel diff = region.endSlices - start;

  // One node per unique voxel corner.
  Eigen::Tensor<NodeId, DIMENSION> nodeIds(diff.x() + 1, diff.y() + 1,
                                           diff.z() + 1);
  std::cout << "GOT DIFF " << nodeIds.dimensions() << std::endl;

  nodeIds.setConstant(NO_NODE);

  // --------------------------------------------------
  // 1. Build vertices and hexahedra.
  // --------------------------------------------------

  const std::array<Voxel, 8> offsets = {
      Voxel(0, 0, 0), Voxel(1, 0, 0), Voxel(1, 1, 0), Voxel(0, 1, 0),
      Voxel(0, 0, 1), Voxel(1, 0, 1), Voxel(1, 1, 1), Voxel(0, 1, 1)};

  for (Voxel v : region) {

    const Voxel relPos = v - start;

    Hex h;
    size_t localNodeId = 0;
    for (const Voxel &offset : offsets) {
      const Voxel vertexPos = relPos + offset;

      NodeId &id = nodeIds(vertexPos.x(), vertexPos.y(), vertexPos.z());

      if (id == NO_NODE) {
        id = static_cast<NodeId>(vertices.size());
        Vertex vertex = vertexPos.cast<double>();
        for (Axis a = 0; a < DIMENSION; a++)
          vertex[a] /= diff[a];
        vertices.push_back(vertex * 2 - Vertex::Ones());
        constrainedFlag.push_back(false);
      }
      assert(id < vertices.size());
      h[localNodeId++] = id;
    }

    assert(localNodeId == 8);

    hexes.push_back(h);
  }

  // --------------------------------------------------
  // 2. Build node adjacency.
  // --------------------------------------------------

  constexpr std::array<std::array<size_t, 2>, 12> edges = {{
      {{0, 1}},
      {{1, 2}},
      {{2, 3}},
      {{3, 0}},
      {{4, 5}},
      {{5, 6}},
      {{6, 7}},
      {{7, 4}},
      {{0, 4}},
      {{1, 5}},
      {{2, 6}},
      {{3, 7}},
  }};
  std::vector<std::unordered_set<NodeId>> neighbours(vertices.size());

  for (const Hex &h : hexes) {
    for (const auto &[a, b] : edges) {
      neighbours[h[a]].insert(h[b]);
      neighbours[h[b]].insert(h[a]);
    }
  }
  // --------------------------------------------------
  // 3. Laplacian smoothing.
  // --------------------------------------------------

  constexpr size_t iterations = 100;

  for (size_t iteration = 0; iteration < iterations; ++iteration) {

    // Apply boundary and external face constraints.
    applyConstraints(nodeIds, start, constrainedFlag, false);

    // Jacobi iteration: read from vertices and write
    // into nextVertices.
    auto nextVertices = vertices;

    for (NodeId id = 0; id < vertices.size(); ++id) {

      if (constrainedFlag[id])
        continue;

      for (int axis = 0; axis < DIMENSION; ++axis) {

        double sum = 0.0;

        for (NodeId neighbourId : neighbours[id]) {
          sum += vertices[neighbourId][axis];
        }

        nextVertices[id][axis] = sum / neighbours[id].size();
      }
    }

    vertices.swap(nextVertices);
  }

  // Reapply constraints after the final smoothing iteration.
  applyConstraints(nodeIds, start, constrainedFlag, true);
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
  std::cout << "VAL NB" << nb << " OUT OF " << constrainedFlag.size()
            << std::endl;

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
