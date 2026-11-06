#include "edge_graph.hpp"
#include <Eigen/Eigen>
#include <cmath>
#include <vector>

void EdgeGraph::buildGraph() {
  adjacency.clear();
  adjacency.resize(mesh.nbBoundaryEdge);
  std::map<std::array<vtkIdType, 2>, size_t> edgeMap;
  size_t edgeId = 0;
  for (size_t faceId = 0; faceId < mesh.faces.size(); faceId++) {
    Face &f = mesh.faces[faceId];
    if (f.adjacentCells[1] != UNDEFINED_ELEMENT)
      continue;
    auto edges = mesh.faceToEdges[faceId];
    for (size_t i = 0; i < NB_NODE_PER_FACE; i++) {
      for (size_t j = 0; j < NB_NODE_PER_FACE; j++) {
        if (i == j)
          continue;
        adjacency[edges[i]].push_back(edges[j]);
        adjacentFace[edges[i]].push_back(faceId);
      }
    }
  }
}

AdjacencyWeight computeEdgeWeights(EdgeGraph &edgeGraph, Vertex axis,
                                   double strictnessFactor) {
  axis.normalize();

  TetrahedralMesh &mesh = edgeGraph.mesh;
  std::vector<std::vector<double>> res(edgeGraph.adjacency.size());

  for (MeshEdgeId nodeId = 0; nodeId < edgeGraph.nbNode(); nodeId++) {
    Vertex p = mesh.edges[nodeId].computeCenter(mesh);
    std::vector<double> nodeWeights;

    for (size_t i = 0; i < edgeGraph.adjacency[nodeId].size(); i++) {
      size_t adjacentNodeId = edgeGraph.adjacency[nodeId][i];
      size_t adjFaceId = edgeGraph.adjacentFace[nodeId][i];
      Vertex padj = mesh.edges[adjacentNodeId].computeCenter(mesh);

      Vertex d = padj - p;
      Vertex faceNormal = mesh.faces[adjFaceId].computeNormal(mesh);
      Vertex v = d.cross(faceNormal).normalized();

      double cosine = v.dot(-axis);

      double w = acos(cosine);
      w = pow(w, strictnessFactor);

      nodeWeights.push_back(w);
    }

    res[nodeId] = nodeWeights;
  }

  return res;
}
