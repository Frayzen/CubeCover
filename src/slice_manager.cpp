#include "slice_manager.hpp"
#include "region.hpp"
#include "voxel_mesh.hpp"
#include <queue>

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
  sliceToConnectedComponents.resize(maxSlice);
  Axis nextAxis = (axis + 1) % 3;
  Axis prevAxis = (axis + 2) % 3;
  size_t nextAxisAmountSlice =
      std::ceil(bboxDiff[nextAxis] / mesh.voxelSize) + 1;
  size_t prevAxisAmountSlice =
      std::ceil(bboxDiff[prevAxis] / mesh.voxelSize) + 1;
  for (SliceId sliceId = 0; sliceId < slices.size(); ++sliceId) {
    Slice &slice = slices[sliceId];

    BoolGrid encounteredFlag(nextAxisAmountSlice, prevAxisAmountSlice);
    encounteredFlag.fill(false);

    for (VoxelId vId : slice) {
      Voxel &v = mesh.voxels[vId];
      encounteredFlag(v[nextAxis], v[prevAxis]) = true;
    }

    BoolGrid visited(nextAxisAmountSlice, prevAxisAmountSlice);
    visited.fill(false);

    for (VoxelId vId : slice) {
      Voxel &v = mesh.voxels[vId];

      size_t x = v[nextAxis];
      size_t y = v[prevAxis];

      if (visited(x, y))
        continue;

      // Found a new connected component.
      sliceToConnectedComponents[sliceId].push_back(ConnectedComponent{vId});

      std::queue<Eigen::Vector2i> queue;
      queue.push({static_cast<int>(x), static_cast<int>(y)});
      visited(x, y) = true;

      while (!queue.empty()) {
        auto p = queue.front();
        queue.pop();

        for (int dx = -1; dx <= 1; ++dx) {
          for (int dy = -1; dy <= 1; ++dy) {
            if (dx == 0 && dy == 0)
              continue;

            int nx = p.x() + dx;
            int ny = p.y() + dy;

            if (nx < 0 || ny < 0 || nx >= nextAxisAmountSlice ||
                ny >= prevAxisAmountSlice)
              continue;

            if (!encounteredFlag(nx, ny) || visited(nx, ny))
              continue;

            visited(nx, ny) = true;
            queue.push({nx, ny});
          }
        }
      }
    }
  }
}

bool SliceManager::isNonConnected(Axis a, Region &r) {
  for (SliceId sliceId = r.startSlices[a]; sliceId < r.endSlices[a];
       sliceId++) {
    size_t amountConnexIn = 0;
    for (ConnectedComponent &cc : sliceToConnectedComponents[sliceId]) {
      if (r.containsVoxel(cc.voxelId))
        amountConnexIn++;
      if (amountConnexIn >= 2)
        break;
    }
    if (amountConnexIn == 1)
      return true;
  }
  return false;
}
