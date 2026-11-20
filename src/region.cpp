#include "region.hpp"
#include "ids.hpp"
#include "region.hpp"
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

std::vector<RegionInterface *>
Region::interfaces(const std::array<bool, DIMENSION> &axisMask) {
  std::vector<RegionInterface *> res;
  for (Axis d = 0; d < DIMENSION; d++) {
    if (!axisMask[d])
      continue;
    for (RegionInterface &interface : children[d])
      res.push_back(&interface);
    // for (RegionInterface &interface : parents[d])
    //   res.push_back(interface);
  }
  return res;
}
