
#include "voxel_mesh.hpp"
#include "constants.hpp"
#include "signed_axis.hpp"

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <vector>

static const VoxelId INVALID_VOXEL_ID = std::numeric_limits<VoxelId>::max();

// -----------------------------------------------------------------------------
// Voxel coordinate helpers
// -----------------------------------------------------------------------------

static Voxel getDimensions(const Eigen::Tensor<VoxelId, DIMENSION> &voxels) {
  return Voxel(static_cast<int>(voxels.dimension(0)),
               static_cast<int>(voxels.dimension(1)),
               static_cast<int>(voxels.dimension(2)));
}

VoxelId VoxelMesh::getId(const Voxel &v) { return voxels(v.x(), v.y(), v.z()); }

Voxel VoxelMesh::getVoxel(VoxelId id) const { return idToVoxel[id]; }

// -----------------------------------------------------------------------------
// Triangle / voxel intersection
// -----------------------------------------------------------------------------

static bool triangleBoxIntersection(const Vertex &A, const Vertex &B,
                                    const Vertex &C, const Vertex &boxMin,
                                    const Vertex &boxMax) {
  const Vertex center = (boxMin + boxMax) * 0.5;
  const Vertex halfSize = (boxMax - boxMin) * 0.5;

  const Vertex v0 = A - center;
  const Vertex v1 = B - center;
  const Vertex v2 = C - center;

  const Vertex edges[3] = {v1 - v0, v2 - v1, v0 - v2};

  const Vertex boxAxes[3] = {Vertex::UnitX(), Vertex::UnitY(), Vertex::UnitZ()};

  auto separated = [&](const Vertex &axis) -> bool {
    if (axis.squaredNorm() < 1e-24)
      return false;

    const double p0 = v0.dot(axis);
    const double p1 = v1.dot(axis);
    const double p2 = v2.dot(axis);

    const double triMin = std::min({p0, p1, p2});
    const double triMax = std::max({p0, p1, p2});

    const double radius = halfSize.dot(axis.cwiseAbs());

    return triMin > radius || triMax < -radius;
  };

  // Box coordinate axes.
  for (const Vertex &axis : boxAxes) {
    if (separated(axis))
      return false;
  }

  // Triangle normal.
  const Vertex normal = edges[0].cross(edges[1]);

  if (separated(normal))
    return false;

  // Cross products of triangle edges and box axes.
  for (const Vertex &edge : edges) {
    for (const Vertex &axis : boxAxes) {
      if (separated(edge.cross(axis)))
        return false;
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
// Voxelize triangle directly into occupancy tensor
// -----------------------------------------------------------------------------

void VoxelMesh::triangleVoxels(const Vertex &A, const Vertex &B,
                               const Vertex &C) {
  if (voxelSize <= 0.0)
    throw std::invalid_argument("voxelSize must be positive");

  const Vertex triMin = A.cwiseMin(B).cwiseMin(C);
  const Vertex triMax = A.cwiseMax(B).cwiseMax(C);

  const Voxel gridMin =
      ((triMin - origin) / voxelSize).array().floor().cast<int>();

  const Voxel gridMax =
      ((triMax - origin) / voxelSize).array().floor().cast<int>();

  const Voxel dimensions = getDimensions(voxels);

  for (int i = std::max(0, gridMin.x());
       i <= std::min(dimensions.x() - 1, gridMax.x()); ++i) {

    for (int j = std::max(0, gridMin.y());
         j <= std::min(dimensions.y() - 1, gridMax.y()); ++j) {

      for (int k = std::max(0, gridMin.z());
           k <= std::min(dimensions.z() - 1, gridMax.z()); ++k) {
        if (voxels(i, j, k) != NO_VOXEL)
          continue;

        const Vertex boxMin = origin + voxelSize * Vertex(i, j, k);

        const Vertex boxMax = boxMin + Vertex::Constant(voxelSize);

        if (triangleBoxIntersection(A, B, C, boxMin, boxMax)) {
          voxels(i, j, k) = nbVoxel++;
          idToVoxel.push_back(Voxel(i, j, k));
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Region normals
// -----------------------------------------------------------------------------

static void updateNormals(Axis axis, Region &oldRegion, Region &newRegion) {
  IntegerSign oldNormal = oldRegion.signedNormals[axis];

  if (oldNormal == UNDEFINED) {
    oldRegion.signedNormals[axis] = NEGATIVE;
    newRegion.signedNormals[axis] = POSITIVE;
  } else if (oldNormal == ZERO) {
    newRegion.signedNormals[axis] = ZERO;
  } else if (oldNormal == NEGATIVE) {
    newRegion.signedNormals[axis] = ZERO;
  } else {
    oldRegion.signedNormals[axis] = ZERO;
  }
}

// -----------------------------------------------------------------------------
// Cut region
// -----------------------------------------------------------------------------

bool VoxelMesh::cut(Axis axis, SliceId splittingSlice) {
  size_t nbRegion = regions.size();
  std::vector<RegionLayerId> newLayerIds(regionGraphs[axis].nbNode());
  std::fill(newLayerIds.begin(), newLayerIds.end(), NO_LAYER);
  for (RegionId rId = 0; rId < nbRegion; rId++) {
    Region &region = regions[rId];
    bool alreadySplittedHere = region.startSlices[axis] == splittingSlice;
    if (alreadySplittedHere)
      return false;
    if (region.containsSlice(axis, splittingSlice)) {
      // Create region
      Region newRegion(region);
      RegionId newRegionId = regions.size();
      newRegion.startSlices[axis] = splittingSlice;
      region.endSlices[axis] = splittingSlice;
      updateNormals(axis, region, newRegion);

      // Layer handling
      RegionLayerId regionLayer = region.layerIds[axis];
      bool encounteredLayer = newLayerIds[regionLayer] != NO_LAYER;
      if (!encounteredLayer) {
        RegionLayerId newLayer = regionGraphs[axis].createNode();
        newLayerIds[regionLayer] = newLayer;
        regionGraphs[axis].adjacency[newLayer] =
            regionGraphs[axis].adjacency[regionLayer];
        regionGraphs[axis].adjacency[regionLayer] = {newLayer};
      }
      RegionLayerId newLayer = newLayerIds[regionLayer];
      newRegion.layerIds[axis] = newLayer;

      // KEEP AT THE END
      regions.push_back(newRegion);
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
// Build voxel mesh
// -----------------------------------------------------------------------------

void VoxelMesh::buildMesh() {

  for (const Vertex &v : tetMesh.vertices)
    bbox.update(v);

  origin = bbox.min;

  const Vertex diff = bbox.max - bbox.min;

  voxelSize = diff.maxCoeff() / 50.0;

  if (voxelSize <= 0.0)
    throw std::runtime_error("Invalid voxel size.");

  const Voxel dimensions = (diff / voxelSize).array().floor().cast<int>() + 1;

  voxels.resize(dimensions.x(), dimensions.y(), dimensions.z());

  voxels.setConstant(NO_VOXEL);

  nbVoxel = 0;
  idToVoxel.clear();
  for (const Face &face : tetMesh.faces) {
    const Vertex &v1 = tetMesh.vertices[face[0]];
    const Vertex &v2 = tetMesh.vertices[face[1]];
    const Vertex &v3 = tetMesh.vertices[face[2]];
    triangleVoxels(v1, v2, v3);
  }

  Region initialRegion(*this);

  initialRegion.signedNormals = {UNDEFINED, UNDEFINED, UNDEFINED};

  initialRegion.layerIds = {0, 0, 0};

  for (Axis d = 0; d < DIMENSION; d++) {
    initialRegion.startSlices[d] = 0;
    initialRegion.endSlices[d] = dimensions[d];
  }
  regions.push_back(initialRegion);
}

// -----------------------------------------------------------------------------
// Find corner
// -----------------------------------------------------------------------------

Voxel VoxelMesh::findCorner(Region &r) {

  Voxel target = Voxel::Zero();

  int numNormals = 0;

  for (Axis d = 0; d < DIMENSION; d++) {

    assert(r.signedNormals[d] != UNDEFINED);

    if (r.signedNormals[d] == ZERO) {
      target[d] = (r.startSlices[d] + r.endSlices[d]) / 2;
      continue;
    }

    numNormals++;

    target[d] =
        r.signedNormals[d] == POSITIVE ? r.endSlices[d] : r.startSlices[d];
  }

  assert(numNormals > 0);

  const Voxel dimensions = getDimensions(voxels);

  double minDist = std::numeric_limits<double>::max();

  Voxel minVoxel = Voxel::Zero();

  for (Voxel v : r) {
    const double dist = (v - target).norm();
    if (dist < minDist) {
      minVoxel = v;
      minDist = dist;
    }
  }

  return minVoxel;
}

// -----------------------------------------------------------------------------
// Shortest path
// -----------------------------------------------------------------------------

std::vector<Voxel> VoxelMesh::findShortestPath(Voxel start, Voxel end, Axis a) {
  using QueueItem = std::pair<double, VoxelId>;

  std::priority_queue<QueueItem, std::vector<QueueItem>,
                      std::greater<QueueItem>>
      queue;

  std::unordered_map<VoxelId, double> distance;
  std::unordered_map<VoxelId, VoxelId> parent;

  const Voxel dimensions = getDimensions(voxels);

  const Vertex A = start.cast<double>();

  const Vertex B = end.cast<double>();

  const Vertex AB = B - A;

  const double lengthSquared = AB.squaredNorm();

  auto distanceToAB = [&](VoxelId id) {
    const Vertex P = getVoxel(id).cast<double>();

    const double t = lengthSquared > 0.0
                         ? std::clamp((P - A).dot(AB) / lengthSquared, 0.0, 1.0)
                         : 0.0;

    const Vertex projection = A + t * AB;

    return (P - projection).norm();
  };

  const std::array<Voxel, 6> directions = {Voxel(1, 0, 0), Voxel(-1, 0, 0),
                                           Voxel(0, 1, 0), Voxel(0, -1, 0),
                                           Voxel(0, 0, 1), Voxel(0, 0, -1)};

  VoxelId startId = getId(start);
  VoxelId endId = getId(end);
  distance[startId] = 0.0;
  parent[startId] = startId;

  queue.push({0.0, startId});

  while (!queue.empty()) {

    auto [currentDistance, current] = queue.top();
    queue.pop();

    if (currentDistance > distance[current])
      continue;

    if (current == endId)
      break;

    const Voxel v = getVoxel(current);

    for (const Voxel &dir : directions) {

      const Voxel next = v + dir;

      if ((next.array() < 0).any() ||
          (next.array() >= dimensions.array()).any()) {
        continue;
      }

      const VoxelId nextId = getId(next);

      if (nextId == INVALID_VOXEL_ID)
        continue;

      const double cost = 1.0 + distanceToAB(nextId);

      const double newDistance = currentDistance + cost;

      if (!distance.contains(nextId) || newDistance < distance[nextId]) {

        distance[nextId] = newDistance;
        parent[nextId] = current;

        queue.push({newDistance, nextId});
      }
    }
  }

  if (!distance.contains(endId))
    throw std::runtime_error("No path found");

  std::vector<Voxel> path;

  for (VoxelId vId = endId; vId != startId; vId = parent[vId])
    path.push_back(getVoxel(vId));

  path.push_back(start);

  std::reverse(path.begin(), path.end());

  return path;
}

// -----------------------------------------------------------------------------
// Find corners
// -----------------------------------------------------------------------------

std::vector<Voxel> VoxelMesh::findCorners() {

  std::vector<Voxel> res;

  for (Region &r : regions)
    res.push_back(findCorner(r));

  return res;
}

// -----------------------------------------------------------------------------
// Find paths between adjacent regions
// -----------------------------------------------------------------------------

std::vector<std::vector<Voxel>> VoxelMesh::findPaths() {

  std::vector<Voxel> corners = findCorners();

  std::vector<std::vector<Voxel>> res;

  const Voxel dimensions = getDimensions(voxels);

  const size_t nbRegion = regions.size();

  for (RegionId rId = 0; rId < nbRegion; rId++) {

    Region &r = regions[rId];

    for (RegionId adjId = rId; adjId < nbRegion; adjId++) {

      Region &adj = regions[adjId];

      Axis axis = 0;
      size_t commonLayers = 0;

      if (adjId == rId)
        continue;

      for (Axis d = 0; d < DIMENSION; d++) {
        if (adj.layerIds[d] == r.layerIds[d])
          commonLayers++;
        else
          axis = d;
      }

      assert(commonLayers < DIMENSION);

      const bool isAdjacent = commonLayers == DIMENSION - 1;

      if (!isAdjacent)
        continue;

      res.push_back(findShortestPath(corners[rId], corners[adjId], axis));
    }
  }

  return res;
}
