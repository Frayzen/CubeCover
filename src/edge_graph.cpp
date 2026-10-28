#include "edge_graph.hpp"
#include <cmath>

void EdgeGraph::buildGraph() {

  std::map<std::array<vtkIdType, 2>, size_t> edgeMap;
  size_t edgeId = 0;
  for (size_t faceId = 0; faceId < mesh.faces.size(); faceId++) {
    Face &f = mesh.faces[faceId];
    size_t edgeIds[3] = {};
    for (size_t i = 0; i < NB_NODE_PER_FACE; i++) {
      std::array<vtkIdType, 2> key = {
          f.vertices[i],
          f.vertices[(i + 1) % NB_NODE_PER_FACE],
      };

      // the key is independent of traversal order
      std::sort(key.begin(), key.end());

      auto it = edgeMap.find(key);
      if (it == edgeMap.end()) {
        edgeMap.emplace(key, edgeId++);
        nodeToVertex.push_back((mesh.vertices[key[0]] + mesh.vertices[key[1]]) /
                               2.0);
        nodeToMeshNodes.push_back(key);
        adjacency.push_back({});
      }
      edgeIds[i] = edgeMap[key];
    }
    for (size_t i = 0; i < NB_NODE_PER_FACE; i++) {
      for (size_t j = 0; j < NB_NODE_PER_FACE; j++) {
        adjacency[edgeIds[i]].push_back({edgeIds[j], faceId});
      }
    }
  }
}

AdjacencyWeight computeWeights(EdgeGraph &graph, Vertex axis,
                               double strictnessFactor) {
  axis.normalize();
  std::vector<std::vector<double>> res;
  for (GraphNodeId nodeId = 0; nodeId < graph.nbNode(); nodeId++) {
    Vertex p = graph.nodeToVertex[nodeId];
    std::vector<double> nodeWeights = {};
    for (auto [adjacentNodeId, adjFaceId] : graph.adjacency[nodeId]) {
      Vertex d = (p - graph.nodeToVertex[adjacentNodeId]) / 2;
      Vertex faceNormal = graph.mesh.faces[adjFaceId].computeNormal(graph.mesh);
      double w = (faceNormal.cross(d)).normalized().dot(axis);
      w = pow(w, strictnessFactor);
      nodeWeights.push_back(w);
    }
    res.push_back(nodeWeights);
  }
  return res;
}
