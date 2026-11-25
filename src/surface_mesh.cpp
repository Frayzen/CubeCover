#include "surface_mesh.hpp"
#include "voxel_mesh.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

SurfaceMesh::SurfaceMesh(const VoxelMesh &mesh)
    : m(mesh), regionToQuadIds(m.regions.size()) {}

// ============================================================
// Find region corner
// ============================================================

NodeId SurfaceMesh::findCorner(DualRegionId regionId) const {

  const DualRegion &r = m.regions[regionId];

  Vertex target = Vertex::Zero();

  int numNormals = 0;

  for (Axis d = 0; d < DIMENSION; ++d) {

    assert(r.boundaryNormals[d] != UNDEFINED);

    if (r.boundaryNormals[d] == ZERO) {

      target[d] = (r.startSlices[d] + r.endSlices[d]) / 2.0;

      continue;
    }

    ++numNormals;

    target[d] =
        r.boundaryNormals[d] == POSITIVE ? r.endSlices[d] : r.startSlices[d];
  }
  target *= m.voxelSize;

  assert(numNormals > 0);

  double minDist = std::numeric_limits<double>::max();
  NodeId minVertex = NO_NODE;

  for (QuadId quadId : regionToQuadIds[regionId]) {
    for (EdgeId edgeId : quadToEdges[quadId]) {
      for (NodeId vertexId : edges[edgeId]) {
        const Vertex &v = nodes[vertexId].cast<double>() * m.voxelSize;
        const double dist = (v - target).norm();

        if (dist < minDist) {
          minDist = dist;
          minVertex = vertexId;
        }
      }
    }
  }

  assert(minVertex != NO_NODE);

  return minVertex;
}

// ============================================================
// Find all region corners
// ============================================================

std::vector<NodeId> SurfaceMesh::findCorners() const {

  std::vector<NodeId> res;

  res.reserve(m.regions.size());

  for (DualRegionId rId = 0; rId < m.regions.size(); ++rId)
    res.push_back(findCorner(rId));

  return res;
}

std::vector<NodeId> SurfaceMesh::getNeighbours(NodeId vertexId) const {

  std::vector<NodeId> res;

  for (EdgeId eId = 0; eId < edges.size(); ++eId) {
    const auto &edge = edges[eId];

    if (edge[0] == vertexId)
      res.push_back(edge[1]);
    else if (edge[1] == vertexId)
      res.push_back(edge[0]);
  }

  return res;
}

// ============================================================
// Dijkstra shortest path
// ============================================================

std::vector<EdgeId> SurfaceMesh::findShortestPath(NodeId start,
                                                  NodeId goal) const {

  using QueueItem = std::pair<double, NodeId>;

  const double INF = std::numeric_limits<double>::max();

  std::vector<double> dist(nodes.size(), INF);

  std::vector<NodeId> previous(nodes.size(), NO_NODE);

  std::vector<EdgeId> previousEdge(nodes.size(), NO_EDGE);

  std::priority_queue<QueueItem, std::vector<QueueItem>,
                      std::greater<QueueItem>>
      queue;

  dist[start] = 0.0;

  queue.emplace(0.0, start);

  while (!queue.empty()) {

    const auto [currentDist, current] = queue.top();

    queue.pop();

    if (currentDist > dist[current])
      continue;

    if (current == goal)
      break;

    for (EdgeId eId : vertexToEdges[current]) {

      const Edge &edge = edges[eId];

      const NodeId neighbor = edge[0] == current ? edge[1] : edge[0];

      assert(current != neighbor);
      Voxel edgeDirection = nodes[current] - nodes[neighbor];
      double weight = edgeDirection.norm();

      std::array<bool, DIMENSION> foundAxis = {false, false, false};
      size_t amountDifferentAxis = 0;
      Axis dirAxis = NO_AXIS;
      for (Axis a = 0; a < DIMENSION; a++)
        if (std::abs(edgeDirection[a]) > m.voxelSize / 2)
          dirAxis = a;
      assert(dirAxis != NO_AXIS);
      for (QuadId qId : edgeToQuads[eId]) {
        Axis a = quadToNormal[qId];
        if (!foundAxis[a])
          amountDifferentAxis++;
        foundAxis[a] = true;
      }

      // Done in order to prefer sharp edges rather than planar surfaces
      if (amountDifferentAxis == 1)
        weight *= (1 + weight / 100.0);

      const double candidate = currentDist + weight;

      if (candidate < dist[neighbor]) {

        dist[neighbor] = candidate;

        previous[neighbor] = current;
        previousEdge[neighbor] = eId;

        queue.emplace(candidate, neighbor);
      }
    }
  }

  if (start != goal && previous[goal] == NO_NODE)
    return {};

  std::vector<EdgeId> path;

  for (NodeId current = goal; current != start; current = previous[current]) {

    if (previousEdge[current] == NO_EDGE)
      return {};

    path.push_back(previousEdge[current]);
  }

  std::reverse(path.begin(), path.end());

  return path;
}

// ============================================================
// Find paths between region corners
// ============================================================

std::vector<SurfacePath> SurfaceMesh::findPaths() {

  regionCorners = findCorners();

  paths.clear();

  for (DualRegionId rId = 0; rId < m.regions.size(); ++rId) {

    const DualRegion &r = m.regions[rId];

    for (Axis a = 0; a < DIMENSION; ++a) {

      for (const DualRegionInterface &interface : r.children[a]) {

        const DualRegionId adjId = interface.regionId;
        const DualRegion &adj = m.regions[adjId];

        SurfacePath path;
        path.regions[0] = rId;
        path.regions[1] = adjId;

        path.normal = interface.normal;

        path.edges = findShortestPath(regionCorners[rId], regionCorners[adjId]);

        const PathId pathId = paths.size();

        for (EdgeId eId : path.edges)
          edgeToPaths[eId].push_back(pathId);

        paths.push_back(std::move(path));
      }
    }
  }

  return paths;
}

// get the square (or box) of adjacent regions. if some adjacent regions are
// missing, sets the resMissingChild flag; otherwise, set it to false
std::set<DualRegionId> SurfaceMesh::getBox(const DualRegionId id,
                                           bool &resMissingChild,
                                           std::set<Axis> excludedAxes) const {
  resMissingChild = false;
  std::set<DualRegionId> dualRegionIds;
  std::function<void(const DualRegion &current, std::set<Axis> excludedAxes)>
      addBox = [&](const DualRegion &current, std::set<Axis> excludedAxes) {
        if (resMissingChild)
          return;

        dualRegionIds.insert(current.id);

        for (Axis axis = 0; axis < DIMENSION; axis++) {
          if (excludedAxes.contains(axis))
            continue;

          if (current.children[axis].empty()) {
            resMissingChild = true;
            return;
          }

          assert(current.children[axis].size() == 1);

          const DualRegion &child =
              m.regions[current.children[axis][0].regionId];

          auto childExcludedAxes = excludedAxes;
          childExcludedAxes.insert(axis);

          addBox(child, std::move(childExcludedAxes));

          if (resMissingChild)
            return;
        }
      };
  addBox(m.regions[id], excludedAxes);
  return dualRegionIds;
}

std::vector<PrimalRegion> SurfaceMesh::buildPrimalRegions() {
  quadToAxisConstraint.clear();
  findPaths();

  if (quadToEdges.size() != quads.size())
    throw std::runtime_error("SurfaceMesh: quadToEdges size mismatch.");

  if (edgeToQuads.size() != edges.size())
    throw std::runtime_error("SurfaceMesh: edgeToQuads size mismatch.");

  if (edgeToPaths.size() != edges.size())
    throw std::runtime_error("SurfaceMesh: edgeToPaths size mismatch.");

  std::vector<bool> visited(quads.size(), false);
  std::vector<Eigen::Vector2i> uv(quads.size());
  std::vector<Axis> colors(quads.size(), NO_AXIS);

  // add boundary surfaces

  for (QuadId seed = 0; seed < quads.size(); ++seed) {

    if (visited[seed])
      continue;

    // Axes of paths bordering this component.
    std::array<size_t, DIMENSION> encounteredAxis{};
    std::set<DualRegionId> encounteredRegions;

    std::vector<QuadId> surfaceQuads;
    std::queue<QuadId> queue;

    visited[seed] = true;
    queue.push(seed);
    DualRegionId seedRegionId = quadToRegion[seed];

    // ----------------------------------------------------
    // Flood fill
    // ----------------------------------------------------

    while (!queue.empty()) {

      const QuadId qId = queue.front();
      queue.pop();

      encounteredRegions.insert(quadToRegion[qId]);

      surfaceQuads.push_back(qId);

      // Inspect the four edges of the current quad.
      for (EdgeId eId : quadToEdges[qId]) {

        if (eId >= edges.size())
          throw std::runtime_error("SurfaceMesh: invalid quad edge index.");

        // ------------------------------------------------
        // Path edge: record its normal and do not cross.
        // ------------------------------------------------

        if (!edgeToPaths[eId].empty()) {

          for (PathId pId : edgeToPaths[eId]) {
            if (pId >= paths.size())
              throw std::runtime_error("SurfaceMesh: invalid path index.");

            if (paths[pId].regions[0] != seedRegionId &&
                paths[pId].regions[1] != seedRegionId)
              continue;

            const Axis a = paths[pId].normal;
            const std::size_t axis = static_cast<std::size_t>(a);

            if (axis >= DIMENSION)
              throw std::runtime_error("SurfaceMesh: invalid path normal.");

            encounteredAxis[axis]++;
          }

          continue;
        }

        // ------------------------------------------------
        // Non-path edge: visit adjacent quads.
        // ------------------------------------------------

        for (QuadId neighbor : edgeToQuads[eId]) {

          if (neighbor >= quads.size())
            throw std::runtime_error(
                "SurfaceMesh: invalid adjacent quad index.");
          if (visited[neighbor])
            continue;

          visited[neighbor] = true;

          queue.push(neighbor);
        }
      }
    }

    assert(encounteredRegions.size() == 4);

    // ----------------------------------------------------
    // Find the missing axis
    // ----------------------------------------------------

    // TODO improve this heuristic
    Axis missingAxis = NO_AXIS;
    size_t minimumCounts = std::numeric_limits<size_t>::max();

    for (Axis a = 0; a < DIMENSION; ++a) {
      if (encounteredAxis[a] >= minimumCounts)
        continue;
      minimumCounts = encounteredAxis[a];
      missingAxis = a;
    }

    // ----------------------------------------------------
    // Assign the missing axis to every quad in the component
    // ----------------------------------------------------

    int normalSignSum = 0;
    for (QuadId qId : surfaceQuads) {
      colors[qId] = missingAxis;
      if (quadToNormal[qId] == missingAxis)
        normalSignSum += quadToNormalSign[qId] == POSITIVE ? 1 : -1;
    }
    SurfaceId surfaceId = surfaces.size();
    IntegerSign averageSign = normalSignSum < 0 ? NEGATIVE : POSITIVE;
    surfaces.push_back({
        .quads = surfaceQuads,
        .normal = missingAxis,
        .sign = averageSign,
        .regions = encounteredRegions,
    });
  }

  // create primal regions

  std::vector<PrimalRegion> primalRegions;
  for (const DualRegion &region : m.regions) {
    bool missingChild = false;
    std::set<DualRegionId> dualRegionIds = getBox(region.id, missingChild);
    if (missingChild)
      continue;
    assert(dualRegionIds.size() == 8);
    PrimalRegion &primalRegion = primalRegions.emplace_back();
    size_t i = 0;
    for (DualRegionId id : dualRegionIds) {
      primalRegion.dualRegions[i] = id;
      primalRegion.corners[i] = nodes[regionCorners[id]];
      i++;
    }
    for (SurfaceId surfaceId = 0; surfaceId < surfaces.size(); surfaceId++) {
      Surface &surface = surfaces[surfaceId];
      if (!surface.allRegionsInSet(dualRegionIds))
        continue;
      SurfaceId &primalSurfaceId =
          primalRegion.getSurfaceId(surface.normal, surface.sign);
      assert(primalSurfaceId == NO_SURFACE);
      primalSurfaceId = surfaceId;
    }
    primalRegion.reorderDualRegionsToVTK(m);
  }

  // add missing internal surfaces
  for (PrimalRegion &primalRegion : primalRegions) {
    for (Axis missingAxis = 0; missingAxis < DIMENSION; missingAxis++) {
      for (IntegerSign sign : {POSITIVE, NEGATIVE}) {
        SurfaceId &primalSurfaceId =
            primalRegion.getSurfaceId(missingAxis, sign);
        if (primalSurfaceId != NO_SURFACE)
          continue;
        std::array<DualRegionId, NB_NODE_PER_QUAD> faceDualRegions =
            primalRegion.getRegionFromFace(missingAxis, sign);

        size_t layer =
            nodes[regionCorners[*faceDualRegions.begin()]][missingAxis];
        for (DualRegionId rId : faceDualRegions)
          assert(nodes[regionCorners[rId]][missingAxis] == layer);
        if (sign == POSITIVE) // because voxel coords are bottom parts
          layer--;
        SurfaceId createdSurfaceId = surfaces.size();
        Surface &createdSurface = surfaces.emplace_back();
        createdSurface.sign = sign;
        createdSurface.normal = missingAxis;
        for (DualRegionId rId : faceDualRegions) {
          Voxel regionNormals = primalRegion.getRegionNormals(rId);
          assert(regionNormals[missingAxis] == (sign == POSITIVE ? 1 : -1));
          createdSurface.regions.insert(rId);
        }
        DualRegionId minRegion = faceDualRegions[0];
        DualRegionId maxRegion = faceDualRegions[2];
        Voxel min = nodes[regionCorners[minRegion]];
        Voxel max = nodes[regionCorners[maxRegion]];
        // TODO fix that
        if (sign == POSITIVE)
          min[missingAxis]--;
        std::cout << "GOT " << min.transpose() << " AND " << max.transpose()
                  << std::endl;
        for (size_t i = min.x(); i < max.x(); i++) {
          for (size_t j = min.y(); j < max.y(); j++) {
            for (size_t k = min.z(); k < max.z(); k++) {
              Voxel v(i, j, k);
              DualRegionId rId = m.voxelToRegion(i, j, k);
              createdSurface.quads.push_back(
                  addQuad(v, rId, missingAxis, sign));
            }
          }
        }

        primalSurfaceId = createdSurfaceId;
      }
    }
  }

  auto dim = m.voxels.dimensions();
  Eigen::Tensor<PrimalRegionId, DIMENSION> voxelToPrimalRegion(dim[0], dim[1],
                                                               dim[2]);
  voxelToPrimalRegion.setConstant(NO_PRIMAL_REGION);
  for (PrimalRegionId prId = 0; prId < primalRegions.size(); prId++) {
    for (SurfaceId sId : primalRegions[prId].surfaceIds) {
      assert(sId != NO_SURFACE);
      Surface &surface = surfaces[sId];
      for (QuadId qId : surface.quads) {
        Voxel v = nodes[quads[qId][0]];
        if (quadToNormalSign[qId] == POSITIVE)
          v[quadToNormal[qId]]--;
        assert(v.x() >= 0 && v.y() >= 0 && v.z() >= 0);
        assert(v.x() < dim[0] && v.y() < dim[1] && v.z() < dim[2]);
        if (voxelToPrimalRegion(v.x(), v.y(), v.z()) != NO_PRIMAL_REGION)
          continue;
        voxelToPrimalRegion(v.x(), v.y(), v.z()) = prId;
        primalRegions[prId].voxels.push_back(v);
      }
    }
  }

  quadToAxisConstraint = colors;
  return primalRegions;
}
