
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

static Voxel getDimensions(const Eigen::Tensor<VoxelId, DIMENSION> &voxels) {
  return Voxel(static_cast<int>(voxels.dimension(0)),
               static_cast<int>(voxels.dimension(1)),
               static_cast<int>(voxels.dimension(2)));
}

VoxelId VoxelMesh::getId(const Voxel &v) { return voxels(v.x(), v.y(), v.z()); }

Voxel VoxelMesh::getVoxel(VoxelId id) const { return idToVoxel[id]; }

static bool tetrahedronBoxIntersection(const Vertex &A, const Vertex &B,
                                       const Vertex &C, const Vertex &D,
                                       const Vertex &boxMin,
                                       const Vertex &boxMax) {

  // Box center and half-size
  const Vertex boxCenter = (boxMin + boxMax) * 0.5;
  const Vertex boxHalfSize = (boxMax - boxMin) * 0.5;

  // Translate tetrahedron vertices relative to box center
  const Vertex v[4] = {A - boxCenter, B - boxCenter, C - boxCenter,
                       D - boxCenter};

  // Tetrahedron edges
  const int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

  constexpr double EPSILON = 1e-12;

  // Test a separating axis
  auto separated = [&](const Vertex &axis) -> bool {
    if (axis.squaredNorm() < EPSILON)
      return false;

    // Project tetrahedron onto axis
    double minProj = v[0].dot(axis);
    double maxProj = minProj;

    for (int i = 1; i < 4; ++i) {
      const double p = v[i].dot(axis);
      minProj = std::min(minProj, p);
      maxProj = std::max(maxProj, p);
    }

    // Project box onto axis
    const double radius = boxHalfSize.x() * std::abs(axis.x()) +
                          boxHalfSize.y() * std::abs(axis.y()) +
                          boxHalfSize.z() * std::abs(axis.z());

    // No intersection if projections do not overlap
    return minProj > radius || maxProj < -radius;
  };

  // 1. Box face normals (X, Y, Z)
  if (separated(Vertex::UnitX()) || separated(Vertex::UnitY()) ||
      separated(Vertex::UnitZ()))
    return false;

  // 2. Tetrahedron face normals
  const Vertex faceNormals[4] = {
      (v[1] - v[0]).cross(v[2] - v[0]), (v[1] - v[0]).cross(v[3] - v[0]),
      (v[2] - v[0]).cross(v[3] - v[0]), (v[2] - v[1]).cross(v[3] - v[1])};

  for (const Vertex &axis : faceNormals) {
    if (separated(axis))
      return false;
  }

  // 3. Cross products between tetrahedron edges
  // and box coordinate axes
  const Vertex boxAxes[3] = {Vertex::UnitX(), Vertex::UnitY(), Vertex::UnitZ()};

  for (const auto &edge : edges) {
    const Vertex e = v[edge[1]] - v[edge[0]];

    for (const Vertex &axis : boxAxes) {
      const Vertex crossAxis = e.cross(axis);

      if (separated(crossAxis))
        return false;
    }
  }

  // No separating axis found: intersection exists
  return true;
}

void VoxelMesh::tetrahedraVoxels(const Vertex &A, const Vertex &B,
                                 const Vertex &C, const Vertex &D) {
  if (voxelSize <= 0.0)
    throw std::invalid_argument("voxelSize must be positive");

  // Tetrahedron bounding box
  const Vertex tetMin = A.cwiseMin(B).cwiseMin(C).cwiseMin(D);

  const Vertex tetMax = A.cwiseMax(B).cwiseMax(C).cwiseMax(D);

  // Convert bounding box to grid coordinates
  const Voxel gridMin =
      ((tetMin - origin) / voxelSize).array().floor().cast<int>();

  const Voxel gridMax =
      ((tetMax - origin) / voxelSize).array().floor().cast<int>();

  const Voxel dimensions = getDimensions(voxels);

  // Iterate over candidate voxels
  for (int i = std::max(0, gridMin.x());
       i <= std::min(dimensions.x() - 1, gridMax.x()); ++i) {

    for (int j = std::max(0, gridMin.y());
         j <= std::min(dimensions.y() - 1, gridMax.y()); ++j) {

      for (int k = std::max(0, gridMin.z());
           k <= std::min(dimensions.z() - 1, gridMax.z()); ++k) {

        if (voxels(i, j, k) != NO_VOXEL)
          continue;

        // Voxel bounding box
        const Vertex boxMin = origin + voxelSize * Vertex(i, j, k);

        const Vertex boxMax = boxMin + Vertex::Constant(voxelSize);

        // Tetrahedron-box intersection test
        if (tetrahedronBoxIntersection(A, B, C, D, boxMin, boxMax)) {
          voxels(i, j, k) = nbVoxel++;
          idToVoxel.push_back(Voxel(i, j, k));
        }
      }
    }
  }
}

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
  for (const Tetrahedron &tet : tetMesh.tets) {
    const Vertex &v1 = tetMesh.vertices[tet[0]];
    const Vertex &v2 = tetMesh.vertices[tet[1]];
    const Vertex &v3 = tetMesh.vertices[tet[2]];
    const Vertex &v4 = tetMesh.vertices[tet[3]];
    tetrahedraVoxels(v1, v2, v3, v4);
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

std::vector<Voxel> VoxelMesh::findCorners() {

  std::vector<Voxel> res;

  for (Region &r : regions)
    res.push_back(findCorner(r));

  return res;
}

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
