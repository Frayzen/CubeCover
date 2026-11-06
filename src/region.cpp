#include "region.hpp"
#include "slice_manager.hpp"

bool Region::containsVoxel(VoxelId voxelId) {
  return std::find(begin(), end(), voxelId) != end();
}

bool Region::containsSlice(Axis axis, SliceId sliceId) {
  return startSlices[axis] <= sliceId && endSlices[axis] > sliceId;
}
