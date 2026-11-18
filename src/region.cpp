#include "region.hpp"
#include "ids.hpp"
#include "region.hpp"
#include "voxel_mesh.hpp"
#include <cmath>

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

void Region::addCut(Region &newRegion, Axis axis,
                    std::set<RegionId> &targetedNeighbours) {

  updateBounds();
  newRegion.updateBounds();

  std::vector<RegionInterface> newRegionChildren;
  std::vector<RegionInterface> regionChildren;

  for (RegionInterface &i : newRegion.children[axis]) {
    Region &child = meshPtr->regions.at(i.regionId);
    assert(child.id == i.regionId);
    if (targetedNeighbours.contains(i.regionId)) {
      regionChildren.push_back(i);
    } else {
      newRegionChildren.push_back(i);
    }
  }

  newRegion.children[axis] = newRegionChildren;
  for (RegionInterface &i : newRegion.children[axis]) {
    Region &child = meshPtr->regions.at(i.regionId);
    for (RegionInterface &i2 : child.parents[axis]) {
      if (i2.regionId != id)
        continue;
      i2.regionId = newRegion.id;
    }
  }
  newRegion.parents[axis] = {RegionInterface(axis, POSITIVE, id)};

  children[axis] = regionChildren;
  children[axis].push_back(RegionInterface(axis, NEGATIVE, newRegion.id));
}

void Region::addCutOld(Region &newRegion, Axis axis, SliceId splittingSlice) {

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
      if (function(*this, child)) {
        Region &next = meshPtr->regions.at(child.regionId);
        next.traverseChildren(function, axisMask, processingFn);
      }
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
      if (function(*this, parent)) {
        Region &next = meshPtr->regions.at(parent.regionId);
        next.traverseChildren(function, axisMask, processingFn);
      }
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
