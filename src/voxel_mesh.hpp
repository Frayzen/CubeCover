#ifndef VOXEL_MESH_HPP
#define VOXEL_MESH_HPP

#include "graph.hpp"
#include "mesh_loader.hpp"
#include "region.hpp"
#include "signed_axis.hpp"
#include "slice_manager.hpp"
#include "vertex.hpp"
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
  VoxelId findCorner(Region &r);
  std::vector<VoxelId>
  findShortestPath(Eigen::Tensor<VoxelId, DIMENSION> &voxelMap, VoxelId start,
                   VoxelId end);
  Graph voxelGraph;

public:
  TetrahedralMesh &tetMesh;
  BBox bbox;

  std::array<SliceManager, DIMENSION> sliceManagers;
  std::vector<Region> regions;

  // A node here is a LayerId
  std::array<Graph, DIMENSION> regionGraphs;

  void buildMesh();
  void splitByConnectedRegions();
  std::vector<Region> getConnectedRegions(RegionId rid);

  double voxelSize;
  Vertex origin;
  std::vector<Voxel> voxels;
  VoxelMesh(TetrahedralMesh &tetMesh)
      : tetMesh(tetMesh), sliceManagers({0, 1, 2}), regionGraphs({1, 1, 1}),
        voxelGraph(0) {
    buildMesh();
  }

  void cut(Axis axis, SliceId slice);
  std::vector<VoxelId> findCorners();
  std::vector<std::vector<VoxelId>> findPaths();

  friend class SliceManager;
};

#endif
