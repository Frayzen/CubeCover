#ifndef SURFACE_MESH_HPP
#define SURFACE_MESH_HPP
#include "ids.hpp"
#include "voxel_mesh.hpp"
#include <vector>

struct SurfacePath {
  std::vector<EdgeId> edges;
  Axis normal;
};

enum FacePathFlag : unsigned char {
  NO_PATH = 0,
  LEFT_PATH = 1 << 0,
  RIGHT_PATH = 1 << 1
};

using Edge = std::array<NodeId, 2>;
using Quad = std::array<NodeId, 4>;
using QuadEdges = std::array<EdgeId, 4>;

struct SurfaceMesh {

  const VoxelMesh &m;

  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  std::vector<Quad> quads;

  // Quad connectivity
  std::vector<QuadEdges> quadToEdges;
  std::vector<RegionId> quadToRegion;
  std::vector<Axis> quadToNormal;
  std::vector<std::vector<QuadId>> regionToQuads;

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

  NodeId findCorner(RegionId regionId) const;

  std::vector<NodeId> findCorners() const;

  std::vector<NodeId> getNeighbours(NodeId vertexId) const;

  std::vector<EdgeId>
  findShortestPath(NodeId start, NodeId goal,
                   std::array<double, DIMENSION> axisPenality) const;

  std::vector<SurfacePath> findPaths();

  void exportMesh(std::string filename) const;
  std::vector<Axis> colorize() const;
};

SurfaceMesh generateSurface(const VoxelMesh &m);
#endif
