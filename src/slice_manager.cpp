#include "slice_manager.hpp"
#include "region.hpp"
#include "voxel_mesh.hpp"

using BoolGrid = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;

void SliceManager::build(VoxelMesh &mesh) {
  auto bboxDiff = (mesh.bbox.max - mesh.bbox.min);
  size_t amountSlice = std::ceil(bboxDiff[axis] / mesh.voxelSize) + 1;

  slices = std::vector<std::vector<VoxelId>>(amountSlice);

  SliceId maxSlice = 0;
  for (VoxelId voxelId = 0; voxelId < mesh.voxels.size(); voxelId++) {
    Voxel &v = mesh.voxels[voxelId];
    SliceId sliceId = v[axis];
    slices[sliceId].push_back(voxelId++);
    maxSlice = std::max(maxSlice, sliceId);
  }
  slices.resize(maxSlice);
}
