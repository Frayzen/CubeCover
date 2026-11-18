#include "region.hpp"
#include "ids.hpp"
#include "region.hpp"
#include "voxel_mesh.hpp"
#include <cmath>
#include <queue>

bool Region::containsSlice(Axis axis, SliceId sliceId) const {
  return startSlices[axis] <= sliceId && endSlices[axis] > sliceId;
}

void Region::updateBounds() {
  Eigen::Vector3<SliceId> newStartSlices;
  newStartSlices.setConstant(MAX_SLICE);
  Eigen::Vector3<SliceId> newEndSlices;
  newEndSlices.setConstant(MIN_SLICE);
  for (Voxel v : *this) {
    newStartSlices = newStartSlices.cwiseMin(v.cast<SliceId>());
    newEndSlices = newEndSlices.cwiseMax(v.cast<SliceId>());
  }
  startSlices = newStartSlices;
  endSlices = newEndSlices + Eigen::Vector3<SliceId>::Ones();
}

void Region::addCut(Region &newRegion, Axis axis, SliceId splittingSlice) {
  updateBounds();
  newRegion.updateBounds();

  newRegion.parents[axis] = {RegionInterface(axis, POSITIVE, id)};
  newRegion.children[axis] = children[axis];

  for (RegionInterface &i : newRegion.children[axis]) {
    Region &child = meshPtr->regions.at(i.regionId);
    for (RegionInterface &i2 : child.parents[axis]) {
      if (i2.regionId != id)
        continue;
      i2.regionId = newRegion.id;
    }
  }

  children[axis] = {RegionInterface(axis, NEGATIVE, newRegion.id)};
}

void Region::traverseChildren(const TraversalFn &function,
                              const std::array<bool, DIMENSION> &axisMask,
                              const ProcessingFn &processingFn) {
  if (processingFn != nullptr && !processingFn(*this))
    return;
  for (Axis d = 0; d < DIMENSION; d++) {
    if (!axisMask[d])
      continue;
    for (RegionInterface &child : children[d]) {
      Region &next = meshPtr->regions.at(child.regionId);
      if (function(*this, child))
        next.traverseChildren(function, axisMask, processingFn);
    }
  }
}
void Region::traverseParents(const TraversalFn &function,
                             const std::array<bool, DIMENSION> &axisMask,
                             const ProcessingFn &processingFn) {
  if (processingFn != nullptr && !processingFn(*this))
    return;
  for (Axis d = 0; d < DIMENSION; d++) {
    if (!axisMask[d])
      continue;
    for (RegionInterface &parent : parents[d]) {
      Region &next = meshPtr->regions.at(parent.regionId);
      if (function(*this, parent))
        next.traverseParents(function, axisMask, processingFn);
    }
  }
}
// TODO change this so that there is no need for traverse
void Region::updateAdjacency(std::vector<RegionId> oldRegionToNewRegion,
                             const std::array<bool, DIMENSION> &axisMask) {
  auto regionFn = [&](Region &r) {
    return oldRegionToNewRegion[r.id] != NO_REGION;
  };
  auto adjFn = [&](Region &r, RegionInterface &interface) {
    if (oldRegionToNewRegion[interface.regionId] != NO_REGION)
      interface.regionId = oldRegionToNewRegion[interface.regionId];
    return true;
  };
  traverseChildren(adjFn, axisMask);
  traverseParents(adjFn, axisMask);
}

std::vector<RegionInterface>
Region::interfaces(const std::array<bool, DIMENSION> &axisMask) {
  std::vector<RegionInterface> res;
  for (Axis d = 0; d < DIMENSION; d++) {
    if (!axisMask[d])
      continue;
    for (RegionInterface &interface : children[d])
      res.push_back(interface);
    for (RegionInterface &interface : parents[d])
      res.push_back(interface);
  }
  return res;
}

Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>
Region::getSliceMatrix(Axis a, SliceId sliceId) {
  assert(containsSlice(a, sliceId));
  Axis nextAxis = (a + 1) % DIMENSION;
  Axis prevAxis = (a + 2) % DIMENSION;

  SliceId startU = startSlices[nextAxis];
  SliceId endU = endSlices[nextAxis];
  size_t diffU = endU - startU;
  SliceId startV = startSlices[prevAxis];
  SliceId endV = endSlices[prevAxis];
  size_t diffV = endV - startV;

  Voxel voxel;
  voxel(a) = sliceId;
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> mat(diffU, diffV);
  for (SliceId u = startU; u < endU; u++) {
    voxel(nextAxis) = u;
    for (SliceId v = startV; v < endV; v++) {
      voxel(prevAxis) = v;
      mat(u - startU, v - startV) = meshPtr->getId(voxel) != NO_VOXEL;
    }
  }

  // Check whether each boundary has a neighbouring voxel
  // outside the region.
  for (SliceId u = startU; u < endU; ++u) {
    voxel(nextAxis) = u;

    voxel(prevAxis) = startV - 1;
    if (startV > 0 && meshPtr->getId(voxel) != NO_VOXEL)
      mat(u - startU, 0) = true;

    voxel(prevAxis) = endV;
    if (endV < meshPtr->voxels.dimension(prevAxis) &&
        meshPtr->getId(voxel) != NO_VOXEL)
      mat(u - startU, diffV - 1) = true;
  }

  for (SliceId v = startV; v < endV; ++v) {
    voxel(prevAxis) = v;

    voxel(nextAxis) = startU - 1;
    if (startU > 0 && meshPtr->getId(voxel) != NO_VOXEL)
      mat(0, v - startV) = true;

    voxel(nextAxis) = endU;
    if (endU < meshPtr->voxels.dimension(nextAxis) &&
        meshPtr->getId(voxel) != NO_VOXEL)
      mat(diffU - 1, v - startV) = true;
  }

  return mat;
}

Eigen::Vector2i Region::getCoordsInSlice(Axis axis, const Voxel &voxel) const {
  assert(axis >= 0 && axis < DIMENSION);

  for (Axis d = 0; d < DIMENSION; ++d)
    assert(containsSlice(d, voxel[d]));

  const Voxel relative = voxel - startSlices;

  const Axis nextAxis = (axis + 1) % DIMENSION;
  const Axis prevAxis = (axis + 2) % DIMENSION;

  return Eigen::Vector2i(relative[nextAxis], relative[prevAxis]);
}

Voxel Region::getVoxelInSlice(Axis axis, const Eigen::Vector2i &coords,
                              SliceId sliceId) const {
  assert(axis >= 0 && axis < DIMENSION);

  Voxel voxel = Voxel::Zero();
  voxel[axis] = sliceId;

  const Axis nextAxis = (axis + 1) % DIMENSION;
  const Axis prevAxis = (axis + 2) % DIMENSION;

  voxel[nextAxis] = coords[0];
  voxel[prevAxis] = coords[1];

  return voxel;
}

using VisitedTensor = Eigen::Tensor<bool, DIMENSION>;
using FlagMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
using FlagTensor = Eigen::Tensor<bool, DIMENSION>;

std::pair<std::vector<Voxel>, std::vector<Voxel>>
Region::splitByFlagMatrix(const FlagMatrix &flagMatrix, const SliceId sliceId,
                          const Axis axis) const {
  assert(axis >= 0 && axis < DIMENSION);

  std::vector<Voxel> inside;
  std::vector<Voxel> outside;

  // endSlices is assumed to be exclusive.
  const Eigen::Vector3i dimensions = endSlices - startSlices;

  assert((dimensions.array() > 0).all());

  FlagTensor occupied(dimensions.x(), dimensions.y(), dimensions.z());
  VisitedTensor visited(dimensions.x(), dimensions.y(), dimensions.z());

  occupied.setConstant(false);
  visited.setConstant(false);

  // Mark voxels belonging to this region.
  for (const Voxel &voxel : voxels) {
    const Voxel local = voxel - startSlices;

    assert((local.array() >= 0).all());
    assert((local.array() < dimensions.array()).all());

    occupied(local.x(), local.y(), local.z()) = true;
  }

  std::queue<Voxel> queue;

  // Flood-fill every connected component intersecting the flagged area.
  for (Eigen::Index x = 0; x < flagMatrix.rows(); ++x) {
    for (Eigen::Index y = 0; y < flagMatrix.cols(); ++y) {
      if (!flagMatrix(x, y))
        continue;

      // Reconstruct the voxel on the splitting slice.
      const Eigen::Vector2i coords(x, y);
      const Voxel startVoxel = getVoxelInSlice(axis, coords, sliceId);

      const Voxel startLocal = startVoxel - startSlices;

      // Check whether the starting voxel belongs to the region.
      if ((startLocal.array() < 0).any() ||
          (startLocal.array() >= dimensions.array()).any())
        continue;

      if (!occupied(startLocal.x(), startLocal.y(), startLocal.z()))
        continue;

      if (visited(startLocal.x(), startLocal.y(), startLocal.z()))
        continue;

      // Start a new connected component.
      queue.push(startVoxel);
      visited(startLocal.x(), startLocal.y(), startLocal.z()) = true;

      while (!queue.empty()) {
        const Voxel voxel = queue.front();
        queue.pop();

        inside.push_back(voxel);

        // Explore the 6-connected neighborhood.
        for (Axis neighborAxis = 0; neighborAxis < DIMENSION; ++neighborAxis) {
          for (int direction : {-1, 1}) {
            Voxel neighbor = voxel;
            neighbor[neighborAxis] += direction;

            // Do not cross the splitting plane.
            if ((voxel[axis] < sliceId && neighbor[axis] >= sliceId) ||
                (voxel[axis] >= sliceId && neighbor[axis] < sliceId))
              continue;

            const Voxel local = neighbor - startSlices;

            // Check region bounds.
            if ((local.array() < 0).any() ||
                (local.array() >= dimensions.array()).any())
              continue;

            const auto lx = local.x();
            const auto ly = local.y();
            const auto lz = local.z();

            if (!occupied(lx, ly, lz) || visited(lx, ly, lz))
              continue;

            visited(lx, ly, lz) = true;
            queue.push(neighbor);
          }
        }
      }
    }
  }

  // Collect all voxels not reached by the flood fills.
  for (const Voxel &voxel : voxels) {
    const Voxel local = voxel - startSlices;

    if (!visited(local.x(), local.y(), local.z()))
      outside.push_back(voxel);
  }

  return {inside, outside};
}
