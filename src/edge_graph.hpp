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
  std::vector<Vertex> nodeToVertex;
  std::vector<std::array<MeshNodeId, 2>> nodeToMeshNodes;

  std::vector<std::vector<std::tuple<GraphNodeId, MeshFaceId>>> adjacency;

  TetrahedralMesh &mesh;

  EdgeGraph(TetrahedralMesh &m) : mesh(m) { buildGraph(); }
  size_t nbNode() { return nodeToVertex.size(); };
};

AdjacencyWeight computeWeights(EdgeGraph &graph, Vertex axis,
                               double strictnessFactor = 10);

#endif
