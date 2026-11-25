#include "dual_region.hpp"
#include "ids.hpp"
#include <cmath>

bool DualRegion::containsSlice(Axis axis, SliceId sliceId) const {
  return startSlices[axis] <= sliceId && endSlices[axis] > sliceId;
}

void DualRegion::updateBounds() {
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

std::vector<DualRegionInterface *>
DualRegion::interfaces(const std::array<bool, DIMENSION> &axisMask) {
  std::vector<DualRegionInterface *> res;
  for (Axis d = 0; d < DIMENSION; d++) {
    if (!axisMask[d])
      continue;
    for (DualRegionInterface &interface : children[d])
      res.push_back(&interface);
    // for (DualRegionInterface &interface : parents[d])
    //   res.push_back(interface);
  }
  return res;
}
