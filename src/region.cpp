#include "region.hpp"
#include "ids.hpp"
#include "region.hpp"
#include "region_iterator.hpp"
#include "voxel_mesh.hpp"
#include <cmath>

bool Region::containsSlice(Axis axis, SliceId sliceId) const {
  return startSlices[axis] <= sliceId && endSlices[axis] > sliceId;
}

RegionIterator Region::begin() const { return RegionIterator(this); }

RegionIterator Region::end() const { return RegionIterator(this, true); }

bool Region::contains(const Voxel &v) const {
  for (Axis d = 0; d < DIMENSION; d++) {
    if (!containsSlice(d, v[d]))
      return false;
  }
  return true;
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
void Region::traverse(const TraversalFn &function,
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
        next.traverse(function, axisMask, processingFn);
    }
    for (RegionInterface &child : children[d]) {
      Region &next = meshPtr->regions.at(child.regionId);
      if (function(*this, child))
        next.traverse(function, axisMask, processingFn);
    }
  }
}
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

Eigen::Vector2i Region::getCoordsInSlice(Axis a, const Voxel &v) {
  assert(contains(v));
  auto relative = v - startSlices;
  Axis nextAxis = (a + 1) % DIMENSION;
  Axis prevAxis = (a + 2) % DIMENSION;
  Eigen::Vector2i res;
  res[0] = relative[nextAxis];
  res[1] = relative[prevAxis];
  return res;
}
