#ifndef VOXEL_MESH_HPP
#define VOXEL_MESH_HPP

#include "mesh_loader.hpp"
#include "region.hpp"
#include "signed_axis.hpp"
#include "slice_manager.hpp"
#include "vertex.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <cmath>
#include <vector>

struct BBox {
  Vertex min = {INFINITY, INFINITY, INFINITY};
  Vertex max = {-INFINITY, -INFINITY, -INFINITY};

  inline void update(const Vertex &v) {
    for (Axis d = 0; d < DIMENSION; d++) {
      min[d] = std::min(min[d], v[d]);
      max[d] = std::max(max[d], v[d]);
    }
  }
};

struct VoxelMesh {
  TetrahedralMesh &tetMesh;
  BBox bbox;

  std::array<SliceManager, DIMENSION> sliceManagers;
  std::vector<Region> regions;
  void buildMesh();

  double voxelSize;
  Vertex origin;
  std::vector<Voxel> voxels;
  VoxelMesh(TetrahedralMesh &tetMesh)
      : tetMesh(tetMesh), sliceManagers({0, 1, 2}) {
    buildMesh();
  }

  void split(Axis axis, SliceId slice);

  friend class SliceManager;
};

#endif
