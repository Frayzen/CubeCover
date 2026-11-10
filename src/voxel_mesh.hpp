#ifndef VOXEL_MESH_HPP
#define VOXEL_MESH_HPP

#include "graph.hpp"
#include "ids.hpp"
#include "mesh_loader.hpp"
#include "region.hpp"
#include "signed_axis.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <cmath>
#include <unsupported/Eigen/CXX11/Tensor>
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
private:
  Voxel findCorner(Region &r);
  std::vector<Voxel> findShortestPath(Voxel start, Voxel end, Axis a);
  Graph voxelGraph;

public:
  TetrahedralMesh &tetMesh;
  BBox bbox;
  size_t nbVoxel;

  Eigen::Tensor<VoxelId, DIMENSION> voxels;
  std::vector<Voxel> idToVoxel;

  std::vector<Region> regions;

  // A node here is a LayerId
  std::array<Graph, DIMENSION> regionGraphs;

  void buildMesh();

  double voxelSize;
  Vertex origin;
  VoxelMesh(TetrahedralMesh &tetMesh)
      : tetMesh(tetMesh), regionGraphs({1, 1, 1}), voxelGraph(0) {
    buildMesh();
  }

  // return true if the cut happened correctly
  // it could return false if a cut already happened on this slide for this axis
  bool cut(Axis axis, SliceId slice);
  std::vector<Voxel> findCorners();
  std::vector<std::vector<Voxel>> findPaths();

  void triangleVoxels(const Vertex &A, const Vertex &B, const Vertex &C);

  VoxelId getId(const Voxel &v);
  Voxel getVoxel(VoxelId id) const;
};

#endif
