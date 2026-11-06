#ifndef SLICE_MANAGER_HPP
#define SLICE_MANAGER_HPP

#include "signed_axis.hpp"
#include "vertex.hpp"
#include <cstddef>
#include <vector>
using SliceId = size_t;

class VoxelMesh;

using Slice = std::vector<VoxelId>;
class SliceManager {
  Axis axis;
  std::vector<Slice> slices;

public:
  void build(VoxelMesh &mesh);
  SliceManager(Axis axis) : axis(axis) {}
  Slice &operator[](SliceId sliceId) { return slices[sliceId]; }
  size_t nbSlice() const { return slices.size(); }
};

#endif
