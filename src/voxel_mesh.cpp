
#include "voxel_mesh.hpp"
#include "constants.hpp"
#include "ids.hpp"
#include "region.hpp"
#include "signed_axis.hpp"
#include "utils.hpp"

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
    oldRegion.signedNormals[axis] = POSITIVE;
    newRegion.signedNormals[axis] = NEGATIVE;
  } else if (oldNormal == ZERO) {
    newRegion.signedNormals[axis] = ZERO;
  } else if (oldNormal == POSITIVE) {
    newRegion.signedNormals[axis] = ZERO;
  } else {
    oldRegion.signedNormals[axis] = ZERO;
  }
}

Eigen::Vector2i VoxelMesh::getSliceSize(Axis axis) const {
  Eigen::Vector2i size;
  int j = 0;

  for (int i = 0; i < DIMENSION; ++i) {
    if (i != axis)
      size[j++] = voxels.dimension(i);
  }

  return size;
}

Eigen::Vector2i VoxelMesh::getCoordsInSlice(Axis axis,
                                            const Voxel &voxel) const {
  Eigen::Vector2i coords;
  int j = 0;

  for (int i = 0; i < DIMENSION; ++i) {
    if (i != axis)
      coords[j++] = voxel[i];
  }

  return coords;
}

std::pair<FlagMatrix, std::vector<RegionId>>
VoxelMesh::computeCutRegions(Axis axis, Voxel origin) {
  std::vector<bool> visitedRegions(regions.size(), false);
  std::vector<RegionId> encounteredRegions;

  FlagMatrix visited(getSliceSize(axis).x(), getSliceSize(axis).y());
  visited.setConstant(false);

  std::queue<Voxel> queue;
  queue.push(origin);

  Eigen::Vector2i originCoords = getCoordsInSlice(axis, origin);
  visited(originCoords.x(), originCoords.y()) = true;

  while (!queue.empty()) {
    Voxel voxel = queue.front();
    queue.pop();

    if (getId(voxel) == NO_VOXEL)
      continue;

    RegionId currentRegionId = voxelToRegion(voxel.x(), voxel.y(), voxel.z());
    if (!visitedRegions[currentRegionId]) {
      visitedRegions[currentRegionId] = true;
      encounteredRegions.push_back(currentRegionId);
    }
    for (Axis a = 0; a < DIMENSION; ++a) {
      if (a == axis)
        continue;

      for (int direction : {-1, 1}) {
        Voxel neighbor = voxel;
        neighbor[a] += direction;

        Eigen::Vector2i coords = getCoordsInSlice(axis, neighbor);

        if (coords.x() < 0 || coords.y() < 0 || coords.x() >= visited.rows() ||
            coords.y() >= visited.cols())
          continue;

        if (visited(coords.x(), coords.y()))
          continue;

        visited(coords.x(), coords.y()) = true;

        if (getId(neighbor) != NO_VOXEL)
          queue.push(neighbor);
      }
    }
  }

  return {visited, encounteredRegions};
}

bool VoxelMesh::cut(Axis axis, Voxel origin) {
  assert(getId(origin) != NO_VOXEL);
  SliceId splittingSlice = origin[axis];
  std::array<bool, DIMENSION> mask = {true, true, true};
  mask[axis] = false;
  std::vector<RegionId> oldRegionToNewRegion(regions.size(), NO_REGION);
  auto [flagMatrix, cutRegionIds] = computeCutRegions(axis, origin);
  std::cout << "SPLITTING " << cutRegionIds.size() << " REGIONS " << std::endl;
  RegionId firstNewRegion = regions.size();
  for (RegionId &cutRegionId : cutRegionIds) {
    RegionId newRegionId = regions.size();
    Region &newRegion = regions.emplace_back(regions[cutRegionId]);
    Region &r = regions[cutRegionId];

    oldRegionToNewRegion[r.id] = newRegionId;
    newRegion.id = newRegionId;

    auto [newVoxels, rVoxels] =
        r.splitByFlagMatrix(flagMatrix, splittingSlice, axis);
    std::cout << "OLD VOXEL AMOUNT " << r.voxels.size() << std::endl;
    std::cout << "INSIDE VOXEL AMOUNT " << newVoxels.size() << std::endl;
    std::cout << "OUTSIDE VOXEL AMOUNT " << rVoxels.size() << std::endl;

    for (Voxel v : newVoxels)
      voxelToRegion(v.x(), v.y(), v.z()) = newRegionId;
    r.voxels = rVoxels;
    newRegion.voxels = newVoxels;
    r.addCut(newRegion, axis, splittingSlice);
  }
  oldRegionToNewRegion.resize(regions.size(), NO_REGION);
  for (RegionId rId = firstNewRegion; rId < regions.size(); rId++) {
    regions[rId].updateAdjacency(oldRegionToNewRegion, mask);
  }
  return true;
}

void VoxelMesh::buildMesh() {
  BBox bbox;
  for (const Vertex &v : tetMesh.vertices)
    bbox.update(v);

  origin = bbox.min;

  const Vertex diff = bbox.max - bbox.min;

  voxelSize = diff.maxCoeff() / 50.0;

  if (voxelSize <= 0.0)
    throw std::runtime_error("Invalid voxel size.");

  const Voxel dimensions = (diff / voxelSize).array().floor().cast<int>() + 1;

  voxels.resize(dimensions.x(), dimensions.y(), dimensions.z());
  voxelToRegion.resize(dimensions.x(), dimensions.y(), dimensions.z());

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

  Region initialRegion(0, *this);

  initialRegion.signedNormals = {UNDEFINED, UNDEFINED, UNDEFINED};
  initialRegion.layerIds = {0, 0, 0};

  for (Axis d = 0; d < DIMENSION; d++) {
    initialRegion.startSlices[d] = 0;
    initialRegion.endSlices[d] = dimensions[d];
  }

  for (size_t i = 0; i < voxels.dimension(0); i++)
    for (size_t j = 0; j < voxels.dimension(1); j++)
      for (size_t k = 0; k < voxels.dimension(2); k++)
        if (voxels(i, j, k) != NO_VOXEL) {
          initialRegion.voxels.push_back(Voxel(i, j, k));
          voxelToRegion(i, j, k) = 0;
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
  auto getNeighbours = [&](Voxel v) {
    std::vector<VoxelId> res;
    for (const Voxel &dir : directions) {
      const Voxel next = v + dir;

      if ((next.array() < 0).any() ||
          (next.array() >= dimensions.array()).any()) {
        continue;
      }

      const VoxelId nextId = getId(next);

      if (nextId == INVALID_VOXEL_ID)
        continue;
      res.push_back(nextId);
    }
    return res;
  };

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

    for (const VoxelId &neighourId : getNeighbours(v)) {
      const Voxel &neighbour = getVoxel(neighourId);
      const double nbNeighbour = getNeighbours(neighbour).size();
      const double cost = distanceToAB(neighourId) + nbNeighbour * nbNeighbour;

      const double newDistance = currentDistance + cost;

      if (!distance.contains(neighourId) ||
          newDistance < distance[neighourId]) {

        distance[neighourId] = newDistance;
        parent[neighourId] = current;

        queue.push({newDistance, neighourId});
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

      bool isAdjacent = false;
      for (Axis d = 0; d < DIMENSION; d++) {
        if (adj.layerIds[d] == r.layerIds[d])
          commonLayers++;
        else
          isAdjacent = adj.startSlices[d] == r.endSlices[d] ||
                       adj.endSlices[d] == r.startSlices[d];
      }
      assert(commonLayers < DIMENSION);
      isAdjacent &= commonLayers == 2;
      if (!isAdjacent)
        continue;

      res.push_back(findShortestPath(corners[rId], corners[adjId], axis));
    }
  }

  return res;
}
