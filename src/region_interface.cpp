#include "region_interface.hpp"
#include "region.hpp"
#include "voxel_mesh.hpp"
void RegionInterface::updateRanges(Region &first, Region &second) {
  SliceId leftFace = NO_SLICE;
  SliceId rightFace = NO_SLICE;
  if (first.startSlices[normal] == second.endSlices[normal]) {
    leftFace = first.startSlices[normal];
    rightFace = first.startSlices[normal] - 1;
  } else if (first.endSlices[normal] == second.startSlices[normal]) {
    leftFace = first.endSlices[normal];
    rightFace = first.endSlices[normal] - 1;
  }

  assert(leftFace != NO_SLICE);
  assert(rightFace != NO_SLICE);

  for (Axis axis = 0; axis < DIMENSION; ++axis) {
    if (axis == normal)
      continue;

    Axis otherAxis = axis == nextAxis ? prevAxis : nextAxis;
    auto &array = axis == nextAxis ? rangesNextAxis : rangesPrevAxis;
    array.clear();

    SliceId startSlice =
        std::max(first.startSlices[axis], second.startSlices[axis]);
    SliceId endSlice = std::min(first.endSlices[axis], second.endSlices[axis]);

    SliceId startV =
        std::max(first.startSlices[otherAxis], second.startSlices[otherAxis]);
    SliceId endV =
        std::min(first.endSlices[otherAxis], second.endSlices[otherAxis]);

    Voxel leftVoxel;
    leftVoxel[normal] = leftFace;

    Voxel rightVoxel;
    rightVoxel[normal] = rightFace;

    bool lastLayerFound = false;
    SliceId beginRangeSlice = NO_SLICE;

    for (SliceId currentSlice = startSlice; currentSlice < endSlice;
         currentSlice++) {
      bool foundForThisLayer = false;
      leftVoxel[axis] = currentSlice;
      rightVoxel[axis] = currentSlice;
      for (SliceId v = startV; v < endV; v++) {
        leftVoxel[otherAxis] = v;
        rightVoxel[otherAxis] = v;
        if (first.meshPtr->getId(leftVoxel) == NO_VOXEL ||
            first.meshPtr->getId(rightVoxel) == NO_VOXEL)
          continue;
        foundForThisLayer = true;
        break;
      }
      if (!foundForThisLayer) {
        if (lastLayerFound)
          array.emplace_back(beginRangeSlice, currentSlice - 1);
      } else if (!lastLayerFound)
        beginRangeSlice = currentSlice;
      lastLayerFound = foundForThisLayer;
    }
    if (lastLayerFound)
      array.emplace_back(beginRangeSlice, endSlice - 1);
  }
}
