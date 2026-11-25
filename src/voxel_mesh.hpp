#ifndef VOXEL_MESH_HPP
#define VOXEL_MESH_HPP

#include "dual_region.hpp"
#include "dual_region_interface.hpp"
#include "graph.hpp"
#include "ids.hpp"
#include "mesh_loader.hpp"
#include "signed_axis.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <cmath>
#include <set>
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

enum class CutFlag : uint8_t { TopRegion, BottomRegion, Undefined };

using CutFlagTensor = Eigen::Tensor<CutFlag, 3>;

struct VoxelPath {
  std::vector<Voxel> nodes;
  Axis normal;
};

struct VoxelMesh {
private:
  Voxel findCorner(DualRegion &r);
  std::vector<Voxel> findShortestPath(Voxel start, Voxel end, Axis a);
  std::vector<SliceId>
  computeLayerCandidates(Axis axis, std::vector<DualRegionId> regions);
  Graph voxelGraph;

  Eigen::Vector2i getSliceSize(Axis axis) const;

  Eigen::Vector2i getCoordsInSlice(Axis axis, const Voxel &voxel) const;

public:
  TetrahedralMesh &tetMesh;
  size_t nbVoxel;

  Eigen::Tensor<VoxelId, DIMENSION> voxels;
  Eigen::Tensor<DualRegionId, DIMENSION> voxelToRegion;
  Eigen::Tensor<bool, DIMENSION> voxelBoundaryFlag;
  Eigen::Tensor<bool, DIMENSION> voxelEdgeFlag;
  Eigen::Tensor<Axis, DIMENSION> voxelPathFlag;
  Eigen::Tensor<std::array<IntegerSign, DIMENSION>, DIMENSION> voxelNormals;
  std::vector<Voxel> idToVoxel;

  std::vector<DualRegion> regions;

  void buildMesh();
  double voxelSize;
  Vertex origin;
  VoxelMesh(TetrahedralMesh &tetMesh) : tetMesh(tetMesh), voxelGraph(0) {
    buildMesh();
  }

  bool isCutAllowed(Axis a, const Voxel &v);
  bool findWay(Axis axis, SliceId sliceId, const DualRegionInterface &first,
               const DualRegionInterface &second);

  // return true if the cut happened correctly
  // it could return false if a cut already happened on this slide for this axis
  bool cut(Axis axis, Voxel origin);
  std::vector<Voxel> findCorners();
  std::vector<VoxelPath> findPaths();

  void tetrahedraVoxels(const Vertex &A, const Vertex &B, const Vertex &C,
                        const Vertex &D);
  Voxel bounds() const {
    Voxel res;
    for (Axis d = 0; d < DIMENSION; d++)
      res[d] = voxels.dimension(d);
    return res;
  };

  void updateBoundaryNormals();

  VoxelId getId(const Voxel &v) const;
  Voxel getVoxel(VoxelId id) const;

  std::tuple<CutFlagTensor, std::set<DualRegionId>, std::set<DualRegionId>,
             std::set<DualRegionId>>
  computeCutRegions(Axis axis, Voxel origin) const;

  std::vector<Voxel> getNeighbours(const Voxel &v,
                                   std::array<bool, DIMENSION> axisMask = {
                                       true, true, true}) const;

  void saveDotfile(Axis a, std::string filename);
};

#endif
