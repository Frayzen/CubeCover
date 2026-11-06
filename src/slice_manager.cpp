#include "slice_manager.hpp"
#include "voxel_mesh.hpp"

void SliceManager::build(VoxelMesh &mesh) {
  double bboxDiff = (mesh.bbox.max - mesh.bbox.min)[axis];
  size_t amountSlice = std::ceil(bboxDiff / mesh.voxelSize) + 1;

  slices = std::vector<std::vector<VoxelId>>(amountSlice);
  sliceConnexParts = std::vector<std::vector<VoxelId>>(amountSlice);

  SliceId maxSlice = 0;

  Axis nextAxis = (axis + 1) % 3;
  Axis prevAxis = (axis + 2) % 3;

  SliceId previousSlice = amountSlice;

  for (VoxelId voxelId = 0; voxelId < mesh.voxels.size(); voxelId++) {
    Voxel &v = mesh.voxels[voxelId];

    SliceId sliceId = v[axis];
    slices[sliceId].push_back(voxelId++);

    bool gapWithPrevious = voxelId != 0;
    if (!gapWithPrevious) {
      Voxel &prev = mesh.voxels[voxelId - 1];
      Voxel dist = prev - v;
      size_t distOtherAxis =
          std::abs(dist[nextAxis]) + std::abs(dist[prevAxis]);
      gapWithPrevious = dist[axis] || distOtherAxis > 1;
    }
    if (gapWithPrevious)
      sliceConnexParts[sliceId].push_back(voxelId);

    maxSlice = std::max(maxSlice, sliceId);
  }
  slices.resize(maxSlice);
}
