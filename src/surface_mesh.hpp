#ifndef SURFACE_MESH_HPP
#define SURFACE_MESH_HPP
#include "constants.hpp"
#include "ids.hpp"
#include "primal_region.hpp"
#include "signed_axis.hpp"
#include "voxel_mesh.hpp"
#include <vector>

using GridVertex = std::array<SliceId, 3>;
using EdgeKey = std::array<NodeId, 2>;

struct SurfacePath {
  std::vector<EdgeId> edges;
  Axis normal;
  DualRegionId regions[2];
};

enum FacePathFlag : unsigned char {
  NO_PATH = 0,
  LEFT_PATH = 1 << 0,
  RIGHT_PATH = 1 << 1
};

using Edge = std::array<NodeId, 2>;
using Quad = std::array<NodeId, 4>;
using QuadEdges = std::array<EdgeId, 4>;

struct Surface {
  std::vector<QuadId> quads;
  Axis normal;
  IntegerSign sign;
  std::set<DualRegionId> regions;

  bool allRegionsInSet(const std::set<DualRegionId> &dualRegionIds) const {
    for (DualRegionId id : regions) {
      if (!dualRegionIds.contains(id))
        return false;
    }
    return true;
  }
};

struct SurfaceMesh {
private:
  // get the square (or box) of adjacent regions. if some adjacent regions are
  // missing, sets the resMissingChild flag; otherwise, leave it as is

  std::set<DualRegionId> getBox(const DualRegionId id, bool &resMissingChild,
                                std::set<Axis> excludedAxes = {}) const;

  NodeId getVertexId(const GridVertex &p);

  EdgeId getEdgeId(NodeId a, NodeId b);

  QuadId addQuad(const std::array<GridVertex, 4> &corners,
                 DualRegionId regionId, Axis normal, IntegerSign normalSign);
  QuadId addQuad(Voxel v, DualRegionId regionId, Axis normal,
                 IntegerSign normalSign);

  std::map<EdgeKey, EdgeId> edgeToId;
  std::map<GridVertex, NodeId> vertexToId;

public:
  const VoxelMesh &m;

  std::vector<Voxel> nodes;
  std::vector<Edge> edges;
  std::vector<Quad> quads;

  // Quad connectivity
  std::vector<QuadEdges> quadToEdges;
  std::vector<DualRegionId> quadToRegion;
  std::vector<Axis> quadToNormal;
  std::vector<IntegerSign> quadToNormalSign;
  std::vector<std::vector<QuadId>> regionToQuadIds;

  std::vector<Axis> quadToAxisConstraint;

  // Primal regions
  std::vector<Surface> surfaces;

  // Edge connectivity
  std::vector<std::vector<QuadId>> edgeToQuads;

  // Vertex connectivity
  std::vector<std::vector<EdgeId>> vertexToEdges;

  // Region corners
  std::vector<NodeId> regionCorners;

  // Paths between region corners
  std::vector<SurfacePath> paths;

  // For each edge, the paths using it.
  std::vector<std::vector<PathId>> edgeToPaths;

  explicit SurfaceMesh(const VoxelMesh &mesh);

  NodeId findCorner(DualRegionId regionId) const;

  std::vector<NodeId> findCorners() const;

  std::vector<NodeId> getNeighbours(NodeId vertexId) const;

  std::vector<EdgeId> findShortestPath(NodeId start, NodeId goal) const;

  std::vector<SurfacePath> findPaths();

  void exportMesh(std::string filename) const;
  std::vector<PrimalRegion> buildPrimalRegions();

  friend SurfaceMesh generateSurface(const VoxelMesh &m);
};

SurfaceMesh generateSurface(const VoxelMesh &m);
#endif
