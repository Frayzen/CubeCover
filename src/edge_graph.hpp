#ifndef EDGE_GRAPH_HPP
#define EDGE_GRAPH_HPP

#include "mesh_loader.hpp"

using MeshNodeId = vtkIdType;
using MeshFaceId = size_t;
using GraphEdgeId = size_t;
using GraphNodeId = size_t;
using AdjacencyWeight = std::vector<std::vector<double>>;

/*
 * Here, each node (of the graph) corresponds to an edge in the mesh
 * Each node is adjacent to another node if their associated edges on the mesh
 * share a face
 */
class EdgeGraph {
private:
  void buildGraph();

public:
  std::vector<std::vector<std::pair<GraphNodeId, MeshFaceId>>> adjacency;

  TetrahedralMesh &mesh;

  EdgeGraph(TetrahedralMesh &m) : mesh(m), adjacency(m.edges.size()) {
    buildGraph();
  }
  size_t nbNode() { return adjacency.size(); };
  double getWeight(AdjacencyWeight &w, size_t from, size_t to);
};

AdjacencyWeight computeWeights(EdgeGraph &graph, Vertex axis,
                               double strictnessFactor = 10);

std::vector<size_t> shortestPath(EdgeGraph &g, AdjacencyWeight &weights,
                                 size_t startingNode);

#endif
