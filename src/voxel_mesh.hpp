#ifndef VOXEL_MESH_HPP
#define VOXEL_MESH_HPP

#include "graph.hpp"
#include "ids.hpp"
#include "mesh_loader.hpp"
#include "region.hpp"
#include "region_interface.hpp"
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

struct VoxelMesh {
private:
  Voxel findCorner(Region &r);
  std::vector<Voxel> findShortestPath(Voxel start, Voxel end, Axis a);
  std::vector<SliceId> computeLayerCandidates(Axis axis,
                                              std::vector<RegionId> regions);
  Graph voxelGraph;

  Eigen::Vector2i getSliceSize(Axis axis) const;

  Eigen::Vector2i getCoordsInSlice(Axis axis, const Voxel &voxel) const;

  std::vector<Voxel> getNeighbours(const Voxel &v,
                                   std::array<bool, DIMENSION> axisMask = {
                                       true, true, true});

public:
  TetrahedralMesh &tetMesh;
  size_t nbVoxel;

  Eigen::Tensor<VoxelId, DIMENSION> voxels;
  Eigen::Tensor<RegionId, DIMENSION> voxelToRegion;
  Eigen::Tensor<bool, DIMENSION> voxelBoundaryFlag;
  Eigen::Tensor<bool, DIMENSION> voxelEdgeFlag;
  std::vector<Voxel> idToVoxel;

  std::vector<Region> regions;

  void buildMesh();
  double voxelSize;
  Vertex origin;
  VoxelMesh(TetrahedralMesh &tetMesh) : tetMesh(tetMesh), voxelGraph(0) {
    buildMesh();
  }

  bool isCutAllowed(Axis a, const Voxel &v);
  bool findWay(Axis axis, SliceId sliceId, const RegionInterface &first,
               const RegionInterface &second);

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

  void updateNormals();

  VoxelId getId(const Voxel &v) const;
  Voxel getVoxel(VoxelId id) const;

  std::tuple<CutFlagTensor, std::set<RegionId>, std::set<RegionId>,
             std::set<RegionId>>
  computeCutRegions(Axis axis, Voxel origin) const;

  void saveDotfile(Axis a, std::string filename);
};

#endif
