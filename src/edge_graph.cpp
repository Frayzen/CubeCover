#include "edge_graph.hpp"
#include <Eigen/Eigen>
#include <cmath>
#include <limits>
#include <queue>
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
        adjacency[edges[i]].push_back(std::make_pair(edges[j], faceId));
      }
    }
  }
}

double EdgeGraph::getWeight(AdjacencyWeight &w, size_t from, size_t to) {
  for (size_t i = 0; i < adjacency[from].size(); i++) {
    if (adjacency[from][i].first == to)
      return w[from][i];
  }
  throw std::invalid_argument("To not in from adjacency list");
}

AdjacencyWeight computeWeights(EdgeGraph &graph, Vertex axis,
                               double strictnessFactor) {
  axis.normalize();

  TetrahedralMesh &mesh = graph.mesh;
  std::vector<std::vector<double>> res(graph.adjacency.size());

  for (MeshEdgeId nodeId = 0; nodeId < graph.nbNode(); nodeId++) {
    Vertex p = mesh.edges[nodeId].computeCenter(mesh);
    std::vector<double> nodeWeights;

    for (auto [adjacentNodeId, adjFaceId] : graph.adjacency[nodeId]) {
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

std::vector<size_t> shortestPath(EdgeGraph &g, AdjacencyWeight &weights,
                                 size_t startingNode) {
  using Weight = double;

  const size_t n = g.nbNode();
  const Weight INF = std::numeric_limits<Weight>::infinity();

  if (startingNode >= n)
    throw std::invalid_argument("Starting node is out of range.");

  if (weights.size() != n)
    throw std::invalid_argument("Weight and graph dimensions do not match.");

  // Validate adjacency and weight dimensions
  for (size_t u = 0; u < n; ++u) {
    if (weights[u].size() != g.adjacency[u].size())
      throw std::invalid_argument(
          "Weight and adjacency dimensions do not match.");

    for (auto [v, faceId] : g.adjacency[u]) {
      if (v >= n)
        throw std::invalid_argument(
            "Adjacency contains an invalid node index.");
    }
  }

  std::vector<Weight> dist(n, INF);
  std::vector<size_t> parent(n, n);

  dist[startingNode] = 0.0;

  // Bellman-Ford: relax all edges at most n-1 times
  for (size_t iter = 0; iter < n - 1; ++iter) {
    bool updated = false;

    for (size_t u = 0; u < n; ++u) {
      if (!std::isfinite(dist[u]))
        continue;

      for (size_t i = 0; i < g.adjacency[u].size(); ++i) {
        size_t v = std::get<0>(g.adjacency[u][i]);
        Weight w = weights[u][i];

        if (!std::isfinite(w))
          throw std::invalid_argument("Edge weight is not finite.");

        Weight candidate = dist[u] + w;

        if (candidate < dist[v]) {
          dist[v] = candidate;
          parent[v] = u;
          updated = true;
        }
      }
    }

    if (!updated)
      break;
  }

  // Build reverse graph to identify nodes that can reach startingNode
  std::vector<std::vector<size_t>> reverseGraph(n);

  for (size_t u = 0; u < n; ++u) {
    for (auto [v, faceId] : g.adjacency[u]) {
      reverseGraph[v].push_back(u);
    }
  }

  std::vector<bool> canReachStart(n, false);
  std::queue<size_t> q;

  canReachStart[startingNode] = true;
  q.push(startingNode);

  while (!q.empty()) {
    size_t u = q.front();
    q.pop();

    for (size_t v : reverseGraph[u]) {
      if (!canReachStart[v]) {
        canReachStart[v] = true;
        q.push(v);
      }
    }
  }

  // Detect reachable negative cycles
  // Only negative cycles that can eventually return to startingNode
  // make the shortest cycle undefined.
  std::vector<bool> affected(n, false);

  for (size_t u = 0; u < n; ++u) {
    if (!std::isfinite(dist[u]))
      continue;

    for (size_t i = 0; i < g.adjacency[u].size(); ++i) {
      size_t v = g.adjacency[u][i].first;
      Weight w = weights[u][i];

      if (dist[u] + w < dist[v]) {
        affected[v] = true;
        q.push(v);
      }
    }
  }

  // Propagate affected nodes
  while (!q.empty()) {
    size_t u = q.front();
    q.pop();

    for (size_t i = 0; i < g.adjacency[u].size(); ++i) {
      size_t v = g.adjacency[u][i].first;

      if (!affected[v]) {
        affected[v] = true;
        q.push(v);
      }
    }
  }

  // A reachable negative cycle that can return to the start
  // means there is no finite shortest cycle.
  for (size_t u = 0; u < n; ++u) {
    if (affected[u] && canReachStart[u])
      throw std::runtime_error("No finite shortest cycle exists: reachable "
                               "negative cycle detected.");
  }

  // Find the cheapest edge returning to startingNode
  Weight bestDistance = INF;
  size_t bestNode = n;

  for (size_t u = 0; u < n; ++u) {
    if (!std::isfinite(dist[u]))
      continue;

    for (size_t i = 0; i < g.adjacency[u].size(); ++i) {
      size_t v = g.adjacency[u][i].first;

      if (v == startingNode) {
        Weight candidate = dist[u] + weights[u][i];

        if (candidate < bestDistance) {
          bestDistance = candidate;
          bestNode = u;
        }
      }
    }
  }

  if (bestNode == n)
    throw std::runtime_error("No cycle returning to the starting node exists.");

  // Reconstruct path
  std::vector<size_t> path;
  size_t current = bestNode;

  while (current != startingNode) {
    if (current >= n || parent[current] == n)
      throw std::runtime_error("Failed to reconstruct the shortest path.");

    path.push_back(current);
    current = parent[current];
  }

  path.push_back(startingNode);

  std::reverse(path.begin(), path.end());

  // Close the cycle
  path.push_back(startingNode);

  return path;
}
