#include "dual_saver.hpp"
#include "region.hpp"

#include <vtkCellData.h>
#include <vtkIdTypeArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <array>
#include <cstring>
#include <stdexcept>
#include <unordered_map>

using Quads = std::vector<std::array<Vertex, 4>>;
using Cut = std::pair<Quads, size_t>;

// Hash for vertices
struct VertexHash {
  size_t operator()(const Vertex &v) const {
    size_t seed = 0;

    for (int i = 0; i < 3; ++i) {
      size_t h = std::hash<double>{}(v[i]);
      seed ^= h + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    return seed;
  }
};

// Equality for vertices
struct VertexEqual {
  bool operator()(const Vertex &a, const Vertex &b) const {
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
  }
};

static void saveCuts(const std::string &filename,
                     const std::vector<Cut> &cuts) {

  auto points = vtkSmartPointer<vtkPoints>::New();

  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // Map vertices to their VTK point IDs
  std::unordered_map<Vertex, vtkIdType, VertexHash, VertexEqual> vertexToId;

  // Cut ID associated with each quad
  auto ids = vtkSmartPointer<vtkIdTypeArray>::New();
  ids->SetName("Id");
  ids->SetNumberOfComponents(1);

  auto getVertexId = [&](const Vertex &v) -> vtkIdType {
    auto it = vertexToId.find(v);

    if (it != vertexToId.end())
      return it->second;

    vtkIdType id = points->InsertNextPoint(v[0], v[1], v[2]);

    vertexToId.emplace(v, id);

    return id;
  };

  // Insert quadrilaterals
  for (const auto &[quads, cutId] : cuts) {

    for (const auto &quad : quads) {

      vtkIdType pointIds[4];

      for (int i = 0; i < 4; ++i)
        pointIds[i] = getVertexId(quad[i]);

      grid->InsertNextCell(VTK_QUAD, 4, pointIds);

      ids->InsertNextValue(static_cast<vtkIdType>(cutId));
    }
  }

  grid->SetPoints(points);
  grid->GetCellData()->AddArray(ids);

  // Write VTU file
  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  writer->SetFileName(filename.c_str());
  writer->SetInputData(grid);

  if (writer->Write() != 1)
    throw std::runtime_error("Failed to write VTU file: " + filename);
}

void saveDual(const std::string &filename, const VoxelMesh &mesh) {
  std::vector<Cut> cuts;

  const Eigen::Vector3i dimensions(mesh.voxels.dimension(0),
                                   mesh.voxels.dimension(1),
                                   mesh.voxels.dimension(2));

  for (const Region &r : mesh.regions) {

    std::array<Quads, DIMENSION> axisQuads;

    for (const Region::Voxel &v : r) {

      const int i = v.x();
      const int j = v.y();
      const int k = v.z();

      if (mesh.voxels(i, j, k) == NO_VOXEL)
        continue;

      const int coords[3] = {i, j, k};

      // Check the two faces perpendicular to each axis
      for (int axis = 0; axis < DIMENSION; ++axis) {

        const int pos = coords[axis];

        // Negative and positive neighboring voxel coordinates
        for (int sign : {-1, 1}) {

          const int neighborPos = pos + sign;

          // Check whether this face is on a region boundary
          const bool onRegionBoundary =
              (sign < 0 && pos == r.startSlices[axis]) ||
              (sign > 0 && pos + 1 == r.endSlices[axis]);

          if (!onRegionBoundary)
            continue;

          // Exclude domain boundaries
          const int facePosition = (sign < 0) ? pos : pos + 1;

          if (facePosition <= 0 || facePosition >= dimensions[axis])
            continue;

          // Check neighboring voxel
          int ni = i;
          int nj = j;
          int nk = k;

          if (axis == 0)
            ni = neighborPos;
          else if (axis == 1)
            nj = neighborPos;
          else
            nk = neighborPos;

          if (mesh.voxels(ni, nj, nk) == NO_VOXEL)
            continue;

          // Construct the four vertices of the face
          const Vertex p = mesh.origin + mesh.voxelSize * Vertex(i, j, k);

          const double s = mesh.voxelSize;

          Quads &quads = axisQuads[axis];

          if (axis == 0) {
            const double x = p.x() + (sign > 0 ? s : 0.0);

            quads.push_back(
                {Vertex(x, p.y(), p.z()), Vertex(x, p.y() + s, p.z()),
                 Vertex(x, p.y() + s, p.z() + s), Vertex(x, p.y(), p.z() + s)});
          } else if (axis == 1) {
            const double y = p.y() + (sign > 0 ? s : 0.0);

            quads.push_back(
                {Vertex(p.x(), y, p.z()), Vertex(p.x(), y, p.z() + s),
                 Vertex(p.x() + s, y, p.z() + s), Vertex(p.x() + s, y, p.z())});
          } else {
            const double z = p.z() + (sign > 0 ? s : 0.0);

            quads.push_back(
                {Vertex(p.x(), p.y(), z), Vertex(p.x() + s, p.y(), z),
                 Vertex(p.x() + s, p.y() + s, z), Vertex(p.x(), p.y() + s, z)});
          }
        }
      }
    }

    // Store one Cut per normal axis
    for (int axis = 0; axis < DIMENSION; ++axis) {
      if (!axisQuads[axis].empty())
        cuts.emplace_back(std::move(axisQuads[axis]), axis);
    }
  }

  saveCuts(filename, cuts);
}
