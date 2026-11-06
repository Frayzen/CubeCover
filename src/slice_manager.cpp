#include "slice_manager.hpp"
#include "voxel_mesh.hpp"

void SliceManager::build(VoxelMesh &mesh) {
  double bboxDiff = (mesh.bbox.max - mesh.bbox.min)[axis];
  size_t amountSlice = std::ceil(bboxDiff / mesh.voxelSize) + 1;
  slices = std::vector<std::vector<VoxelId>>(amountSlice);
  VoxelId id = 0;
  SliceId maxSlice = 0;
  for (Voxel &v : mesh.voxels) {
    SliceId cur = v[axis];
    slices[cur].push_back(id++);
    maxSlice = std::max(maxSlice, cur);
  }
  slices.resize(maxSlice);
}
