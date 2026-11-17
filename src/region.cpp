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

bool Region::isCutAllowed(Axis axis, SliceId splittingSlice) {
  for (Axis d = 0; d < DIMENSION; d++) {
    if (axis == d)
      continue;
    bool canCut = false;
    for (RegionInterface &interface : children[d]) {
      if (interface.contains(axis, splittingSlice))
        canCut = true;
    }
    for (RegionInterface &interface : parents[d])
      if (interface.contains(axis, splittingSlice))
        canCut = true;
    if (!canCut)
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
  newRegion.parents[axis][0].updateRanges(*this, newRegion);
  newRegion.children[axis] = children[axis];

  for (RegionInterface &i : newRegion.children[axis]) {
    Region &child = meshPtr->regions.at(i.regionId);
    i.updateRanges(child, newRegion);
    for (RegionInterface &i2 : child.parents[axis]) {
      if (i2.regionId != id)
        continue;
      i2.regionId = newRegion.id;
      i2.updateRanges(newRegion, child);
    }
  }

  children[axis] = {RegionInterface(axis, NEGATIVE, newRegion.id)};
  children[axis][0].updateRanges(*this, newRegion);
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
