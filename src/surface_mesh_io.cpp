#include "surface_mesh.hpp"
#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkLine.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertex.h>
#include <vtkXMLUnstructuredGridWriter.h>

// ============================================================
// Generate surface mesh
// ============================================================

SurfaceMesh generateSurface(const VoxelMesh &m) {

  SurfaceMesh surface(m);

  using GridVertex = std::array<SliceId, 3>;
  using EdgeKey = std::array<NodeId, 2>;

  std::map<GridVertex, NodeId> vertexToId;
  std::map<EdgeKey, EdgeId> edgeToId;

  const SliceId nx = m.voxels.dimension(0);
  const SliceId ny = m.voxels.dimension(1);
  const SliceId nz = m.voxels.dimension(2);

  // --------------------------------------------------------
  // Get or create vertex
  // --------------------------------------------------------

  auto getVertexId = [&](const GridVertex &p) -> NodeId {
    auto [it, inserted] = vertexToId.emplace(p, surface.nodes.size());

    if (inserted) {

      Voxel &v = surface.nodes.emplace_back(p[0], p[1], p[2]);
      surface.vertexToEdges.emplace_back();
    }

    return it->second;
  };

  // --------------------------------------------------------
  // Get or create edge
  // --------------------------------------------------------

  auto getEdgeId = [&](NodeId a, NodeId b) -> EdgeId {
    EdgeKey key = {std::min(a, b), std::max(a, b)};

    auto [it, inserted] = edgeToId.emplace(key, surface.edges.size());

    if (inserted) {

      const EdgeId eId = it->second;

      surface.edges.push_back(key);

      surface.edgeToQuads.emplace_back();
      surface.edgeToPaths.emplace_back();

      surface.vertexToEdges[key[0]].push_back(eId);
      surface.vertexToEdges[key[1]].push_back(eId);
    }

    return it->second;
  };

  // --------------------------------------------------------
  // Add quad
  // --------------------------------------------------------

  auto addQuad = [&](const std::array<GridVertex, 4> &corners,
                     RegionId regionId, Axis normal) {
    Quad quad;

    QuadEdges quadEdges;

    for (size_t i = 0; i < 4; ++i)
      quad[i] = getVertexId(corners[i]);

    for (size_t i = 0; i < 4; ++i) {

      const size_t next = (i + 1) % 4;

      quadEdges[i] = getEdgeId(quad[i], quad[next]);
    }

    const QuadId qId = surface.quads.size();

    surface.quads.push_back(quad);
    surface.quadToEdges.push_back(quadEdges);
    surface.quadToRegion.push_back(regionId);
    surface.quadToNormal.push_back(normal);
    surface.regionToQuads[regionId].push_back(qId);

    for (EdgeId eId : quadEdges)
      surface.edgeToQuads[eId].push_back(qId);
  };

  // --------------------------------------------------------
  // Occupancy
  // --------------------------------------------------------

  auto isOccupied = [&](SliceId x, SliceId y, SliceId z) -> bool {
    if (x >= nx || y >= ny || z >= nz)
      return false;

    return m.voxels(x, y, z) != NO_VOXEL;
  };

  // --------------------------------------------------------
  // Generate exposed voxel faces
  // --------------------------------------------------------

  for (SliceId x = 0; x < nx; ++x) {

    for (SliceId y = 0; y < ny; ++y) {

      for (SliceId z = 0; z < nz; ++z) {

        if (!isOccupied(x, y, z))
          continue;

        const RegionId rId = m.voxelToRegion(x, y, z);

        const SliceId xp = x + 1;
        const SliceId yp = y + 1;
        const SliceId zp = z + 1;

        // LEFT (-X)
        if (x == 0 || !isOccupied(x - 1, y, z)) {

          addQuad({{{{x, y, z}}, {{x, y, zp}}, {{x, yp, zp}}, {{x, yp, z}}}},
                  rId, 0);
        }

        // RIGHT (+X)
        if (!isOccupied(xp, y, z)) {

          addQuad(
              {{{{xp, y, z}}, {{xp, yp, z}}, {{xp, yp, zp}}, {{xp, y, zp}}}},
              rId, 0);
        }

        // FRONT (-Y)
        if (y == 0 || !isOccupied(x, y - 1, z)) {

          addQuad({{{{x, y, z}}, {{xp, y, z}}, {{xp, y, zp}}, {{x, y, zp}}}},
                  rId, 1);
        }

        // BACK (+Y)
        if (!isOccupied(x, yp, z)) {

          addQuad(
              {{{{x, yp, z}}, {{x, yp, zp}}, {{xp, yp, zp}}, {{xp, yp, z}}}},
              rId, 1);
        }

        // BOTTOM (-Z)
        if (z == 0 || !isOccupied(x, y, z - 1)) {

          addQuad({{{{x, y, z}}, {{x, yp, z}}, {{xp, yp, z}}, {{xp, y, z}}}},
                  rId, 2);
        }

        // TOP (+Z)
        if (!isOccupied(x, y, zp)) {

          addQuad(
              {{{{x, y, zp}}, {{xp, y, zp}}, {{xp, yp, zp}}, {{x, yp, zp}}}},
              rId, 2);
        }
      }
    }
  }

  return surface;
}

void SurfaceMesh::exportMesh(std::string filename) const {

  if (quads.size() != quadToRegion.size()) {
    throw std::runtime_error("SurfaceMesh: quadToRegion size mismatch.");
  }

  if (regionCorners.size() != m.regions.size()) {
    throw std::runtime_error("SurfaceMesh: regionCorners size mismatch.");
  }

  if (paths.size() == 0 && !m.regions.empty()) {
    // Empty paths are allowed.
  }

  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(static_cast<vtkIdType>(nodes.size()));

  for (NodeId i = 0; i < nodes.size(); ++i) {
    const Vertex &v = nodes[i].cast<double>() * m.voxelSize;

    points->SetPoint(static_cast<vtkIdType>(i), v.x(), v.y(), v.z());
  }

  vtkNew<vtkUnstructuredGrid> grid;
  grid->SetPoints(points);

  vtkNew<vtkIntArray> cellTypes;
  cellTypes->SetName("CellType");

  vtkNew<vtkIntArray> quadColor;
  quadColor->SetName("QuadColor");

  vtkNew<vtkIntArray> regionIds;
  regionIds->SetName("RegionId");

  vtkNew<vtkIntArray> pathIds;
  pathIds->SetName("PathId");

  vtkNew<vtkIntArray> pathNormals;
  pathNormals->SetName("PathNormal");

  // ========================================================
  // Quads
  // ========================================================

  vtkNew<vtkQuad> quad;

#define COLORIZE
  for (QuadId qId = 0; qId < quads.size(); ++qId) {

    const auto &face = quads[qId];

    for (int i = 0; i < NB_NODE_PER_QUAD; ++i) {

      if (face[i] >= nodes.size()) {
        throw std::runtime_error("SurfaceMesh: invalid quad vertex index.");
      }

      quad->GetPointIds()->SetId(i, static_cast<vtkIdType>(face[i]));
    }

    grid->InsertNextCell(VTK_QUAD, quad->GetPointIds());

    cellTypes->InsertNextValue(0);
    regionIds->InsertNextValue(static_cast<int>(quadToRegion[qId]));
#ifdef COLORIZE
    quadColor->InsertNextValue(quadToAxisConstraint[qId]);
#endif

    pathIds->InsertNextValue(-1);
    pathNormals->InsertNextValue(-1);
  }

  // ========================================================
  // Region corners
  // ========================================================

  for (RegionId rId = 0; rId < regionCorners.size(); ++rId) {

    const NodeId vertexId = regionCorners[rId];

    if (vertexId >= nodes.size()) {
      throw std::runtime_error("SurfaceMesh: invalid corner vertex index.");
    }

    vtkNew<vtkVertex> vertex;
    vertex->GetPointIds()->SetId(0, static_cast<vtkIdType>(vertexId));

    grid->InsertNextCell(VTK_VERTEX, vertex->GetPointIds());

    cellTypes->InsertNextValue(1);
    regionIds->InsertNextValue(static_cast<int>(rId));
#ifdef COLORIZE
    quadColor->InsertNextValue(-1);
#endif

    pathIds->InsertNextValue(-1);
    pathNormals->InsertNextValue(-1);
  }

  vtkNew<vtkIntArray> isCorner;
  isCorner->SetName("IsCorner");
  isCorner->SetNumberOfComponents(1);
  isCorner->SetNumberOfTuples(static_cast<vtkIdType>(nodes.size()));

  // Initialize all nodes as non-corners.
  for (vtkIdType i = 0; i < static_cast<vtkIdType>(nodes.size()); ++i)
    isCorner->SetValue(i, 0);

  // Mark region corners.
  for (NodeId vertexId : regionCorners) {
    if (vertexId >= nodes.size())
      throw std::runtime_error("SurfaceMesh: invalid corner vertex index.");

    isCorner->SetValue(static_cast<vtkIdType>(vertexId), 1);
  }

  // ========================================================
  // Paths
  // ========================================================

  for (PathId pId = 0; pId < paths.size(); ++pId) {

    const SurfacePath &path = paths[pId];

    for (EdgeId eId : path.edges) {

      if (eId >= edges.size()) {
        throw std::runtime_error("SurfaceMesh: invalid path edge index.");
      }

      const auto &edge = edges[eId];

      if (edge[0] >= nodes.size() || edge[1] >= nodes.size()) {
        throw std::runtime_error("SurfaceMesh: invalid edge vertex index.");
      }

      vtkNew<vtkLine> line;

      line->GetPointIds()->SetId(0, static_cast<vtkIdType>(edge[0]));

      line->GetPointIds()->SetId(1, static_cast<vtkIdType>(edge[1]));

      grid->InsertNextCell(VTK_LINE, line->GetPointIds());

      cellTypes->InsertNextValue(2);
      regionIds->InsertNextValue(-1);
#ifdef COLORIZE
      quadColor->InsertNextValue(-1);
#endif

      pathIds->InsertNextValue(static_cast<int>(pId));
      pathNormals->InsertNextValue(static_cast<int>(path.normal));
    }
  }

  // ========================================================
  // Cell data
  // ========================================================

  grid->GetCellData()->AddArray(cellTypes);
  grid->GetCellData()->AddArray(regionIds);
  grid->GetCellData()->AddArray(pathIds);
  grid->GetCellData()->AddArray(pathNormals);
#ifdef COLORIZE
  grid->GetCellData()->AddArray(quadColor);
#endif
  grid->GetPointData()->AddArray(isCorner);

  // ========================================================
  // Write VTU
  // ========================================================

  vtkNew<vtkXMLUnstructuredGridWriter> writer;

  writer->SetFileName(filename.c_str());
  writer->SetInputData(grid);
  writer->SetDataModeToBinary();

  if (writer->Write() != 1) {
    throw std::runtime_error("SurfaceMesh: failed to write VTU file: " +
                             filename);
  }
}
