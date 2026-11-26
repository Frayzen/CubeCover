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

// --------------------------------------------------------
// Get or create vertex
// --------------------------------------------------------

NodeId SurfaceMesh::getVertexId(const GridVertex &p) {
  auto [it, inserted] = vertexToId.emplace(p, nodes.size());

  if (inserted) {
    nodes.emplace_back(p[0], p[1], p[2]);
    vertexToEdges.emplace_back();
  }

  return it->second;
}

// --------------------------------------------------------
// Get or create edge
// --------------------------------------------------------

EdgeId SurfaceMesh::getEdgeId(NodeId a, NodeId b) {
  EdgeKey key = {std::min(a, b), std::max(a, b)};

  auto [it, inserted] = edgeToId.emplace(key, edges.size());

  if (inserted) {
    const EdgeId eId = it->second;

    edges.push_back(key);

    edgeToQuads.emplace_back();
    edgeToPaths.emplace_back();

    vertexToEdges[key[0]].push_back(eId);
    vertexToEdges[key[1]].push_back(eId);
  }

  return it->second;
}

// --------------------------------------------------------
// Add quad
// --------------------------------------------------------

QuadId SurfaceMesh::addQuad(const std::array<GridVertex, 4> &corners,
                            DualRegionId regionId, Axis normal,
                            IntegerSign normalSign) {
  Quad quad;
  QuadEdges quadEdges;

  for (size_t i = 0; i < 4; ++i)
    quad[i] = getVertexId(corners[i]);

  for (size_t i = 0; i < 4; ++i) {
    const size_t next = (i + 1) % 4;
    quadEdges[i] = getEdgeId(quad[i], quad[next]);
  }

  const QuadId qId = quads.size();

  quads.push_back(quad);
  quadToEdges.push_back(quadEdges);
  quadToRegion.push_back(regionId);
  quadToNormal.push_back(normal);
  quadToNormalSign.push_back(normalSign);
  regionToQuadIds[regionId].push_back(qId);

  for (EdgeId eId : quadEdges)
    edgeToQuads[eId].push_back(qId);

  return qId;
}

QuadId SurfaceMesh::addQuad(Voxel v, DualRegionId regionId, Axis normal,
                            IntegerSign normalSign) {
  std::array<GridVertex, 4> corners;
  GridVertex cur = {v.x(), v.y(), v.z()};
  if (normalSign == POSITIVE)
    cur[normal]++;
  Axis nextAxis = (normal + 1) % DIMENSION;
  Axis prevAxis = (normal + 2) % DIMENSION;
  size_t cornerId = 0;
  std::array<Eigen::Vector2i, NB_NODE_PER_QUAD> offsets = {
      Eigen::Vector2i(0, 0),
      Eigen::Vector2i(0, 1),
      Eigen::Vector2i(1, 1),
      Eigen::Vector2i(1, 0),
  };
  for (auto offset : offsets) {
    cur[nextAxis] = v[nextAxis] + offset[0];
    cur[prevAxis] = v[prevAxis] + offset[1];
    corners[cornerId++] = cur;
  }
  return addQuad(corners, regionId, normal, normalSign);
}

// ============================================================
// Generate surface mesh
// ============================================================

SurfaceMesh generateSurface(const VoxelMesh &m) {

  SurfaceMesh surface(m);

  const SliceId nx = m.voxels.dimension(0);
  const SliceId ny = m.voxels.dimension(1);
  const SliceId nz = m.voxels.dimension(2);

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

        const DualRegionId rId = m.voxelToRegion(x, y, z);

        const SliceId xp = x + 1;
        const SliceId yp = y + 1;
        const SliceId zp = z + 1;

        // LEFT (-X)
        if (x == 0 || !isOccupied(x - 1, y, z)) {

          surface.addQuad(
              {{{{x, y, z}}, {{x, y, zp}}, {{x, yp, zp}}, {{x, yp, z}}}}, rId,
              0, NEGATIVE);
        }

        // RIGHT (+X)
        if (!isOccupied(xp, y, z)) {

          surface.addQuad(
              {{{{xp, y, z}}, {{xp, yp, z}}, {{xp, yp, zp}}, {{xp, y, zp}}}},
              rId, 0, POSITIVE);
        }

        // FRONT (-Y)
        if (y == 0 || !isOccupied(x, y - 1, z)) {

          surface.addQuad(
              {{{{x, y, z}}, {{xp, y, z}}, {{xp, y, zp}}, {{x, y, zp}}}}, rId,
              1, NEGATIVE);
        }

        // BACK (+Y)
        if (!isOccupied(x, yp, z)) {

          surface.addQuad(
              {{{{x, yp, z}}, {{x, yp, zp}}, {{xp, yp, zp}}, {{xp, yp, z}}}},
              rId, 1, POSITIVE);
        }

        // BOTTOM (-Z)
        if (z == 0 || !isOccupied(x, y, z - 1)) {

          surface.addQuad(
              {{{{x, y, z}}, {{x, yp, z}}, {{xp, yp, z}}, {{xp, y, z}}}}, rId,
              2, NEGATIVE);
        }

        // TOP (+Z)
        if (!isOccupied(x, y, zp)) {

          surface.addQuad(
              {{{{x, y, zp}}, {{xp, y, zp}}, {{xp, yp, zp}}, {{x, yp, zp}}}},
              rId, 2, POSITIVE);
        }
      }
    }
  }

  return surface;
}

void SurfaceMesh::exportMesh(
    std::string filename,
    const std::vector<PrimalRegion> &primalRegions) const {

  if (quads.size() != quadToRegion.size()) {
    throw std::runtime_error("SurfaceMesh: quadToRegion size mismatch.");
  }

  if (regionCorners.size() != m.regions.size()) {
    throw std::runtime_error("SurfaceMesh: regionCorners size mismatch.");
  }

  if (paths.size() == 0 && !m.regions.empty()) {
    // Empty paths are allowed.
  }

  // Map each boundary quad to its primal region and surface.
  std::vector<int> quadToPrimalRegionId(quads.size(), -1);
  std::vector<int> quadToSurfaceId(quads.size(), -1);
  std::vector<int> quadToSurfaceNormal(quads.size(), -1);
  std::vector<int> quadToSurfaceSign(quads.size(), -1);

  for (PrimalRegionId prId = 0; prId < primalRegions.size(); ++prId) {

    for (SurfaceId sId : primalRegions[prId].surfaceIds) {

      if (sId == NO_SURFACE)
        continue;

      if (sId >= surfaces.size())
        throw std::runtime_error("SurfaceMesh: invalid surface index.");

      const Surface &surface = surfaces[sId];

      for (QuadId qId : surface.quads) {

        if (qId >= quads.size())
          throw std::runtime_error("SurfaceMesh: invalid surface quad index.");

        // A quad must have a unique primal-region/surface association.
        if (quadToPrimalRegionId[qId] != -1) {
          if (quadToPrimalRegionId[qId] != static_cast<int>(prId) ||
              quadToSurfaceId[qId] != static_cast<int>(sId)) {
            throw std::runtime_error("SurfaceMesh: quad belongs to multiple "
                                     "primal region surfaces.");
          }
          continue;
        }

        quadToPrimalRegionId[qId] = static_cast<int>(prId);
        quadToSurfaceId[qId] = static_cast<int>(sId);
        quadToSurfaceNormal[qId] = static_cast<int>(surface.normal);
        quadToSurfaceSign[qId] = surface.sign == POSITIVE ? 1 : -1;
      }
    }
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

  vtkNew<vtkIntArray> primalRegionIds;
  primalRegionIds->SetName("PrimalRegionId");

  vtkNew<vtkIntArray> surfaceIds;
  surfaceIds->SetName("SurfaceId");

  vtkNew<vtkIntArray> surfaceNormals;
  surfaceNormals->SetName("SurfaceNormal");

  vtkNew<vtkIntArray> surfaceSigns;
  surfaceSigns->SetName("SurfaceSign");

  vtkNew<vtkIntArray> isBoundary;
  isBoundary->SetName("IsBoundary");

  // ========================================================
  // Quads
  // ========================================================

  vtkNew<vtkQuad> quad;

  if (quadToAxisConstraint.empty())
    std::cerr << "Warning: quad to axis constraint is empty" << std::endl;

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
    if (quadToAxisConstraint.size() > qId)
      quadColor->InsertNextValue(quadToAxisConstraint[qId]);
    else
      quadColor->InsertNextValue(-1);

    pathIds->InsertNextValue(-1);
    pathNormals->InsertNextValue(-1);
    primalRegionIds->InsertNextValue(quadToPrimalRegionId[qId]);
    surfaceIds->InsertNextValue(quadToSurfaceId[qId]);
    surfaceNormals->InsertNextValue(quadToSurfaceNormal[qId]);
    surfaceSigns->InsertNextValue(quadToSurfaceSign[qId]);

    isBoundary->InsertNextValue(quadToPrimalRegionId[qId] != -1 ? 1 : 0);
  }

  // ========================================================
  // Region corners
  // ========================================================

  for (DualRegionId rId = 0; rId < regionCorners.size(); ++rId) {

    const NodeId vertexId = regionCorners[rId];

    if (vertexId >= nodes.size()) {
      throw std::runtime_error("SurfaceMesh: invalid corner vertex index.");
    }

    vtkNew<vtkVertex> vertex;
    vertex->GetPointIds()->SetId(0, static_cast<vtkIdType>(vertexId));

    grid->InsertNextCell(VTK_VERTEX, vertex->GetPointIds());

    cellTypes->InsertNextValue(1);
    regionIds->InsertNextValue(static_cast<int>(rId));
    quadColor->InsertNextValue(-1);

    pathIds->InsertNextValue(-1);
    pathNormals->InsertNextValue(-1);
    primalRegionIds->InsertNextValue(-1);
    surfaceIds->InsertNextValue(-1);
    surfaceNormals->InsertNextValue(-1);
    surfaceSigns->InsertNextValue(-1);
    isBoundary->InsertNextValue(0);
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
      quadColor->InsertNextValue(-1);

      pathIds->InsertNextValue(static_cast<int>(pId));
      pathNormals->InsertNextValue(static_cast<int>(path.normal));
      primalRegionIds->InsertNextValue(-1);
      surfaceIds->InsertNextValue(-1);
      surfaceNormals->InsertNextValue(-1);
      surfaceSigns->InsertNextValue(-1);
      isBoundary->InsertNextValue(0);
    }
  }

  // ========================================================
  // Cell data
  // ========================================================

  grid->GetCellData()->AddArray(cellTypes);
  grid->GetCellData()->AddArray(regionIds);
  grid->GetCellData()->AddArray(pathIds);
  grid->GetCellData()->AddArray(pathNormals);
  grid->GetCellData()->AddArray(quadColor);
  grid->GetPointData()->AddArray(isCorner);
  grid->GetCellData()->AddArray(primalRegionIds);
  grid->GetCellData()->AddArray(surfaceIds);
  grid->GetCellData()->AddArray(surfaceNormals);
  grid->GetCellData()->AddArray(surfaceSigns);
  grid->GetCellData()->AddArray(isBoundary);

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
