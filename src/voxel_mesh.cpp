
#include "voxel_mesh.hpp"
#include "constants.hpp"
#include "ids.hpp"
#include "region.hpp"
#include "signed_axis.hpp"
#include "utils.hpp"
#include <fstream>

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

VoxelId VoxelMesh::getId(const Voxel &v) const {
  return voxels(v.x(), v.y(), v.z());
}

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

  const Voxel dimensions = bounds();

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

void VoxelMesh::updateNormals() {
  std::vector<std::array<bool, DIMENSION>> hasParents(regions.size(),
                                                      {false, false, false});
  for (Region &region : regions) {
    for (Axis a = 0; a < DIMENSION; a++)
      for (RegionInterface &i : region.children[a])
        hasParents[i.regionId][a] = true;
  }
  for (Region &region : regions)
    for (Axis a = 0; a < DIMENSION; a++) {
      if (region.children[a].size() == 0)
        region.signedNormals[a] = NEGATIVE;
      else if (!hasParents[region.id][a])
        region.signedNormals[a] = POSITIVE;
      else
        region.signedNormals[a] = ZERO;
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

std::tuple<CutFlagTensor, std::set<RegionId>, std::set<RegionId>,
           std::set<RegionId>>
VoxelMesh::computeCutRegions(Axis axis, Voxel origin) const {
  assert(axis >= 0 && axis < DIMENSION);

  const Voxel dimensions = bounds();
  const SliceId splittingSlice = origin[axis];

  CutFlagTensor flags(dimensions[0], dimensions[1], dimensions[2]);
  flags.setConstant(CutFlag::Undefined);

  std::set<RegionId> encounteredRegions;
  std::set<RegionId> topRegionNeighbors;
  std::set<RegionId> bottomRegionNeighbors;

  // --------------------------------------------------
  // 1. Find the connected cut using a 2D flood fill.
  // --------------------------------------------------

  FlagMatrix cut(getSliceSize(axis).x(), getSliceSize(axis).y());
  cut.setConstant(false);

  std::queue<Voxel> queue;

  if (voxels(origin.x(), origin.y(), origin.z()) == NO_VOXEL)
    return {flags, encounteredRegions, topRegionNeighbors,
            bottomRegionNeighbors};

  queue.push(origin);

  const Eigen::Vector2i originCoords = getCoordsInSlice(axis, origin);
  cut(originCoords.x(), originCoords.y()) = true;

  while (!queue.empty()) {
    const Voxel voxel = queue.front();
    queue.pop();

    for (Axis neighborAxis = 0; neighborAxis < DIMENSION; ++neighborAxis) {
      if (neighborAxis == axis)
        continue;

      for (int direction : {-1, 1}) {
        Voxel neighbor = voxel;
        neighbor[neighborAxis] += direction;

        const Eigen::Vector2i coords = getCoordsInSlice(axis, neighbor);

        if (coords.x() < 0 || coords.y() < 0 || coords.x() >= cut.rows() ||
            coords.y() >= cut.cols())
          continue;

        if (cut(coords.x(), coords.y()))
          continue;

        if (voxels(neighbor.x(), neighbor.y(), neighbor.z()) == NO_VOXEL)
          continue;

        cut(coords.x(), coords.y()) = true;
        queue.push(neighbor);
      }
    }
  }

  // --------------------------------------------------
  // 2. Flood fill one side of the cut.
  // --------------------------------------------------

  auto floodFill = [&](Voxel seed, CutFlag sideFlag,
                       std::set<RegionId> &sideNeighbors) {
    if ((seed.array() < 0).any() || (seed.array() >= dimensions.array()).any())
      return;

    if (voxels(seed.x(), seed.y(), seed.z()) == NO_VOXEL)
      return;

    // Do not overwrite voxels already assigned to the
    // other side. A valid separating cut should prevent
    // the two traversals from overlapping.
    if (flags(seed.x(), seed.y(), seed.z()) != CutFlag::Undefined)
      return;

    Eigen::Tensor<bool, 3> visited(dimensions[0], dimensions[1], dimensions[2]);
    visited.setConstant(false);

    std::queue<Voxel> sideQueue;

    sideQueue.push(seed);
    visited(seed.x(), seed.y(), seed.z()) = true;
    flags(seed.x(), seed.y(), seed.z()) = sideFlag;

    const RegionId seedRegion = voxelToRegion(seed.x(), seed.y(), seed.z());

    encounteredRegions.insert(seedRegion);

    while (!sideQueue.empty()) {
      const Voxel voxel = sideQueue.front();
      sideQueue.pop();

      const RegionId currentRegionId =
          voxelToRegion(voxel.x(), voxel.y(), voxel.z());

      for (Axis neighborAxis = 0; neighborAxis < DIMENSION; ++neighborAxis) {

        for (int direction : {-1, 1}) {
          Voxel neighbor = voxel;
          neighbor[neighborAxis] += direction;

          if ((neighbor.array() < 0).any() ||
              (neighbor.array() >= dimensions.array()).any())
            continue;

          const auto x = neighbor.x();
          const auto y = neighbor.y();
          const auto z = neighbor.z();

          if (visited(x, y, z))
            continue;

          if (voxels(x, y, z) == NO_VOXEL)
            continue;

          // Do not enter voxels already assigned to
          // the opposite side.
          if (flags(x, y, z) != CutFlag::Undefined)
            continue;

          const RegionId neighborRegionId = voxelToRegion(x, y, z);

          if (neighborRegionId != currentRegionId) {

            // Do not change region along the splitting axis.
            if (neighborAxis == axis) {
              sideNeighbors.insert(neighborRegionId);
              continue;
            }

            // Do not change region outside the projected cut.
            const Eigen::Vector2i coords = getCoordsInSlice(axis, voxel);

            if (!cut(coords.x(), coords.y())) {
              sideNeighbors.insert(neighborRegionId);
              continue;
            }
          }

          // Do not cross the cut.
          if (neighborAxis == axis) {
            const bool crossingCut = (voxel[axis] < splittingSlice &&
                                      neighbor[axis] >= splittingSlice) ||
                                     (voxel[axis] >= splittingSlice &&
                                      neighbor[axis] < splittingSlice);

            if (crossingCut) {
              const Eigen::Vector2i coords = getCoordsInSlice(axis, voxel);

              if (cut(coords.x(), coords.y()))
                continue;
            }
          }

          visited(x, y, z) = true;
          flags(x, y, z) = sideFlag;

          if (neighborRegionId != NO_REGION)
            encounteredRegions.insert(neighborRegionId);

          sideQueue.push(neighbor);
        }
      }
    }
  };

  // --------------------------------------------------
  // 3. Traverse the top and bottom sides independently.
  // --------------------------------------------------

  // Top side: origin and positive axis direction.
  floodFill(origin, CutFlag::TopRegion, topRegionNeighbors);

  // Bottom side: immediately below the splitting slice.
  Voxel bottomOrigin = origin;
  bottomOrigin[axis]--;

  floodFill(bottomOrigin, CutFlag::BottomRegion, bottomRegionNeighbors);

  return {flags, encounteredRegions, topRegionNeighbors, bottomRegionNeighbors};
}

bool VoxelMesh::cut(Axis axis, Voxel origin) {
  assert(getId(origin) != NO_VOXEL);
  SliceId splittingSlice = origin[axis];
  std::array<bool, DIMENSION> mask = {true, true, true};
  mask[axis] = false;
  auto [flags, encounteredRegions, topRegionNeighbors, bottomRegionNeighbors] =
      computeCutRegions(axis, origin);

  // for (auto r : topRegionNeighbors) {
  //   std::cout << "TOP VOLUME NEIGH " << r << std::endl;
  // }

  // for (auto r : bottomRegionNeighbors) {
  //   std::cout << "BOT VOLUME NEIGH " << r << std::endl;
  // }

  std::vector<RegionId> cutRegionToNewRegion(regions.size(), NO_REGION);
  for (RegionId cutRegionId : encounteredRegions) {
    Region &newRegion = regions.emplace_back(regions.size(), *this);
    Region &cutRegion = regions[cutRegionId];
    cutRegionToNewRegion[cutRegionId] = newRegion.id;
    std::vector<Voxel> cutRegionVoxels;
    std::vector<Voxel> newRegionVoxels;
    for (Voxel v : cutRegion) {
      if (flags(v.x(), v.y(), v.z()) == CutFlag::TopRegion) {
        newRegionVoxels.push_back(v);
        voxelToRegion(v.x(), v.y(), v.z()) = newRegion.id;
      } else {
        assert(flags(v.x(), v.y(), v.z()) == CutFlag::BottomRegion);
        cutRegionVoxels.push_back(v);
      }
    }
    cutRegion.voxels = cutRegionVoxels;
    newRegion.voxels = newRegionVoxels;
    cutRegion.updateBounds();
    newRegion.updateBounds();
  }

  for (RegionId cutRegionId : encounteredRegions) {
    Region &newRegion = regions[cutRegionToNewRegion[cutRegionId]];
    Region &cutRegion = regions[cutRegionId];

    newRegion.children[axis] = {RegionInterface(axis, NEGATIVE, cutRegion.id)};
    for (Axis a = 0; a < DIMENSION; a++) {
      std::vector<RegionInterface> newCutRegionChildren;
      for (RegionInterface &i : cutRegion.children[a]) {
        bool partOfCutRegions = encounteredRegions.contains(i.regionId);
        bool partOfTopRegions = partOfCutRegions;
        partOfTopRegions |= topRegionNeighbors.contains(i.regionId);
        bool partOfBottomRegions = partOfCutRegions;
        partOfBottomRegions |= bottomRegionNeighbors.contains(i.regionId);
        if (partOfTopRegions) {
          RegionInterface &created = newRegion.children[a].emplace_back(i);
          if (cutRegionToNewRegion[created.regionId] != NO_REGION)
            created.regionId = cutRegionToNewRegion[i.regionId];
        }
        if (partOfBottomRegions)
          newCutRegionChildren.push_back(i);
      }
      cutRegion.children[a] = newCutRegionChildren;
    }
  }

  for (RegionId neighborId : topRegionNeighbors) {
    for (RegionInterface *i : regions[neighborId].interfaces()) {
      if (cutRegionToNewRegion[i->regionId] != NO_REGION)
        i->regionId = cutRegionToNewRegion[i->regionId];
    }
  }
  return true;
}

std::vector<Voxel>
VoxelMesh::getNeighbours(const Voxel &v, std::array<bool, DIMENSION> axisMask) {
  const Voxel dimensions = bounds();
  std::vector<Voxel> res;
  for (Axis a = 0; a < DIMENSION; a++) {
    if (!axisMask[a])
      continue;
    for (int da : {-1, 1}) {
      Voxel neighbor = v;
      neighbor(a) += da;
      if ((neighbor.array() < 0).any() ||
          (neighbor.array() >= dimensions.array()).any())
        continue;
      if (getId(neighbor) != NO_VOXEL)
        res.push_back(neighbor);
    }
  }
  return res;
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
  voxelBoundaryFlag.resize(dimensions.x(), dimensions.y(), dimensions.z());
  voxelEdgeFlag.resize(dimensions.x(), dimensions.y(), dimensions.z());

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
          voxelBoundaryFlag(i, j, k) =
              getNeighbours(Voxel(i, j, k)).size() != 6;
          voxelEdgeFlag(i, j, k) = getNeighbours(Voxel(i, j, k)).size() <= 4;
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

  const Voxel dimensions = bounds();

  double minDist = std::numeric_limits<double>::max();

  Voxel minVoxel = Voxel::Zero();

  for (Voxel v : r) {
    if (!voxelEdgeFlag(v.x(), v.y(), v.z()))
      continue;
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

  const Voxel dimensions = bounds();

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
    for (const Voxel &neighbor : getNeighbours(v)) {
      if (!voxelBoundaryFlag(neighbor.x(), neighbor.y(), neighbor.z()))
        continue;
      const VoxelId neighborId = getId(neighbor);
      const double nbNeighbour = getNeighbours(neighbor).size();
      const double cost = distanceToAB(neighborId) + nbNeighbour * nbNeighbour;

      const double newDistance = currentDistance + cost;

      if (!distance.contains(neighborId) ||
          newDistance < distance[neighborId]) {

        distance[neighborId] = newDistance;
        parent[neighborId] = current;

        queue.push({newDistance, neighborId});
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

  const Voxel dimensions = bounds();
  const size_t nbRegion = regions.size();

  for (RegionId rId = 0; rId < nbRegion; rId++) {

    Region &r = regions[rId];

    for (Axis a = 0; a < DIMENSION; a++)
      for (RegionInterface interface : r.children[a]) {
        RegionId adjId = interface.regionId;
        Region &adj = regions[adjId];
        res.push_back(
            findShortestPath(corners[rId], corners[adjId], interface.normal));
      }
  }
  return res;
}

bool VoxelMesh::isCutAllowed(Axis a, const Voxel &v) {
  Region &r = regions[voxelToRegion(v.x(), v.y(), v.z())];
  for (Axis d = 0; d < DIMENSION; d++) {
    if (d == a)
      continue;
    std::array<bool, DIMENSION> axisMask = {false, false, false};
    axisMask[d] = true;
    if (r.interfaces(axisMask).size() == 0)
      return false;
  }
  return true;
}

void VoxelMesh::saveDotfile(Axis a, std::string filename) {
  assert(a >= 0 && a < DIMENSION);

  std::ofstream file(filename);

  if (!file.is_open()) {
    throw std::runtime_error("Could not open DOT file: " + filename);
  }

  file << "digraph RegionGraph {\n";
  file << "  rankdir=TB;\n";
  file << "  node [shape=box, style=rounded];\n\n";

  // Write nodes.
  for (const Region &region : regions) {
    file << "  " << region.id << " [label=\"Region " << region.id << "\"];\n";
  }

  file << "\n";

  // Write edges using child relationships.
  for (const Region &region : regions) {
    for (const RegionInterface &child : region.children[a]) {
      file << "  " << region.id << " -> " << child.regionId << ";\n";

      // bool foundParent = false;
      // for (const RegionInterface &parent :
      // regions[child.regionId].parents[a]) {
      //   if (parent.regionId == region.id) {
      //     foundParent = true;
      //     break;
      //   }
      // }
      // assert(foundParent);
    }
  }

  file << "}\n";
}
