#include "region.hpp"
#include "ids.hpp"
#include "voxel_mesh.hpp"
#include <cmath>

bool Region::containsSlice(Axis axis, SliceId sliceId) const {
  return startSlices[axis] <= sliceId && endSlices[axis] > sliceId;
}

bool Region::contains(const Voxel &v) const {
  for (Axis d = 0; d < DIMENSION; d++) {
    if (!containsSlice(d, v[d]))
      return false;
  }
  return true;
}

Region::RegionIterator::RegionIterator(const Region *region, bool atEnd)
    : region_(region), atEnd_(atEnd) {
  if (!region_ || atEnd_) {
    atEnd_ = true;
    return;
  }

  x_ = static_cast<int>(region_->startSlices.x());
  y_ = static_cast<int>(region_->startSlices.y());
  z_ = static_cast<int>(region_->startSlices.z());

  advanceToOccupied();
}

Region::Voxel Region::RegionIterator::operator*() const {
  return Voxel(x_, y_, z_);
}

void Region::RegionIterator::advanceToOccupied() {
  if (!region_ || !region_->meshPtr) {
    atEnd_ = true;
    return;
  }

  const auto &voxels = region_->meshPtr->voxels;

  const int endX = static_cast<int>(region_->endSlices.x());
  const int endY = static_cast<int>(region_->endSlices.y());
  const int endZ = static_cast<int>(region_->endSlices.z());

  while (z_ < endZ) {

    if (x_ >= endX) {
      x_ = static_cast<int>(region_->startSlices.x());
      ++y_;
    }

    if (y_ >= endY) {
      y_ = static_cast<int>(region_->startSlices.y());
      ++z_;
      continue;
    }

    if (x_ >= 0 && y_ >= 0 && z_ >= 0 && x_ < voxels.dimension(0) &&
        y_ < voxels.dimension(1) && z_ < voxels.dimension(2) &&
        voxels(x_, y_, z_) != NO_VOXEL) {
      atEnd_ = false;
      return;
    }

    ++x_;
  }

  atEnd_ = true;
}

Region::RegionIterator &Region::RegionIterator::operator++() {
  if (!atEnd_) {
    ++x_;
    advanceToOccupied();
  }
  return *this;
}

bool Region::RegionIterator::operator==(const RegionIterator &other) const {
  if (atEnd_ && other.atEnd_)
    return region_ == other.region_;

  return region_ == other.region_ && atEnd_ == other.atEnd_ && x_ == other.x_ &&
         y_ == other.y_ && z_ == other.z_;
}

void Region::addCut(Axis axis, SliceId splittingSlice) {
  Eigen::Vector3<SliceId> newStartSlices;
  newStartSlices.setConstant(MAX_SLICE);
  Eigen::Vector3<SliceId> newEndSlices;
  newEndSlices.setConstant(MIN_SLICE);
  size_t i = 0;
  for (Voxel v : *this) {
    newStartSlices = newStartSlices.cwiseMin(v.cast<SliceId>());
    newEndSlices = newEndSlices.cwiseMax(v.cast<SliceId>());
    i++;
  }
  startSlices = newStartSlices;
  endSlices = newEndSlices + Eigen::Vector3<SliceId>::Ones();
}
