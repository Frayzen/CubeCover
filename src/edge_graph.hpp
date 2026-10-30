#ifndef EDGE_GRAPH_HPP
#define EDGE_GRAPH_HPP

#include "graph.hpp"
#include "mesh_loader.hpp"

using MeshNodeId = vtkIdType;
using MeshFaceId = size_t;
using GraphEdgeId = size_t;
using MeshEdgeId = size_t;

/*
 * Here, each node (of the graph) corresponds to an edge in the mesh
 * Each node is adjacent to another node if their associated edges on the mesh
 * share a face
 */
class EdgeGraph : public Graph {
private:
  void buildGraph();

public:
  TetrahedralMesh &mesh;

  std::vector<std::vector<size_t>>
      adjacentFace; // associated with Graph.adjacency

  EdgeGraph(TetrahedralMesh &m)
      : mesh(m), Graph(m.edges.size()), adjacentFace(m.edges.size()) {
    buildGraph();
  }
};

AdjacencyWeight computeEdgeWeights(EdgeGraph &edgeGraph, Vertex axis,
                                   double strictnessFactor = 10);

#endif
