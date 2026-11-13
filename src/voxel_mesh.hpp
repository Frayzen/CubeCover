#ifndef VOXEL_MESH_HPP
#define VOXEL_MESH_HPP

#include "axis_layer_manager.hpp"
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
  std::vector<SliceId> computeLayerCandidates(Axis axis,
                                              std::vector<RegionId> regions);
  Graph voxelGraph;

public:
  TetrahedralMesh &tetMesh;
  size_t nbVoxel;

  Eigen::Tensor<VoxelId, DIMENSION> voxels;
  std::vector<Voxel> idToVoxel;

  std::array<AxisLayerManager, DIMENSION> layerManagers = {
      AxisLayerManager(*this, 0),
      AxisLayerManager(*this, 1),
      AxisLayerManager(*this, 2),
  };
  std::vector<Region> regions;

  void buildMesh();

  double voxelSize;
  Vertex origin;
  VoxelMesh(TetrahedralMesh &tetMesh) : tetMesh(tetMesh), voxelGraph(0) {
    buildMesh();
  }

  // return true if the cut happened correctly
  // it could return false if a cut already happened on this slide for this axis
  bool cut(Axis axis, Voxel origin);
  std::vector<Voxel> findCorners();
  std::vector<std::vector<Voxel>> findPaths();

  void tetrahedraVoxels(const Vertex &A, const Vertex &B, const Vertex &C,
                        const Vertex &D);
  Voxel bounds() const {
    Voxel res;
    for (Axis d = 0; d < DIMENSION; d++)
      res[d] = voxels.dimension(d);
    return res;
  };

  VoxelId getId(const Voxel &v);
  Voxel getVoxel(VoxelId id) const;

  std::vector<RegionId> getAdjacentRegions(RegionId regionId, Axis axis);
};

#endif
