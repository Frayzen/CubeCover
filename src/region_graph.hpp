#ifndef AXIS_GRAPH_HPP
#define AXIS_GRAPH_HPP

#include "edge_graph.hpp"
#include "many_to_one.hpp"

using NodeId = size_t;
using FaceId = size_t;
using EdgeId = size_t;
using RegionId = size_t;
using Axis = size_t;
using Path = std::pair<std::vector<EdgeId>, Axis>;

class RegionGraph {
private:
  EdgeGraph &g;
  TetrahedralMesh &mesh;
  size_t nbRegion;
  ManyToOne<NodeId, RegionId> regions;

  std::vector<Vertex> nodeNormals;
  void computeNodeNormals();

  void fillConnectedFaces(size_t startingFace, std::vector<size_t> edgeIds,
                          size_t axis);

public:
  std::vector<Path> paths;
  RegionGraph(EdgeGraph &g);

  void insertLoop(std::vector<size_t> edgesId, size_t axis);
  size_t getFaceRegion(size_t nodeId) { return regions[nodeId]; };

  std::vector<NodeId> getCandidates(RegionId regionId);

  void saveMesh(std::string filename);
};

#endif
