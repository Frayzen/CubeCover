#include "graph.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <queue>
#include <stdexcept>

double Graph::getWeight(AdjacencyWeight &w, size_t from, size_t to) {
  for (size_t i = 0; i < adjacency[from].size(); i++) {
    if (adjacency[from][i] == to)
      return w[from][i];
  }
  throw std::invalid_argument("To not in from adjacency list");
}

std::vector<size_t> bellmanFordShortestPath(Graph &g, AdjacencyWeight &weights,
                                            size_t startingNode,
                                            size_t endingNode) {
  using Weight = double;

  const size_t n = g.nbNode();
  const Weight INF = std::numeric_limits<Weight>::infinity();

  if (startingNode >= n)
    throw std::invalid_argument("Starting node is out of range.");

  if (endingNode >= n)
    throw std::invalid_argument("Ending node is out of range.");

  if (weights.size() != n)
    throw std::invalid_argument("Weight and graph dimensions do not match.");

  // Validate adjacency and weight dimensions
  for (size_t u = 0; u < n; ++u) {
    if (weights[u].size() != g.adjacency[u].size())
      throw std::invalid_argument(
          "Weight and adjacency dimensions do not match.");

    for (auto v : g.adjacency[u]) {
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
        size_t v = g.adjacency[u][i];
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

  // Build reverse graph to identify nodes that can reach endingNode
  std::vector<std::vector<size_t>> reverseGraph(n);

  for (size_t u = 0; u < n; ++u) {
    for (auto v : g.adjacency[u]) {
      reverseGraph[v].push_back(u);
    }
  }

  std::vector<bool> canReachEnd(n, false);
  std::queue<size_t> q;

  canReachEnd[endingNode] = true;
  q.push(endingNode);

  while (!q.empty()) {
    size_t u = q.front();
    q.pop();

    for (size_t v : reverseGraph[u]) {
      if (!canReachEnd[v]) {
        canReachEnd[v] = true;
        q.push(v);
      }
    }
  }

  // Detect reachable negative cycles
  // Only negative cycles that can eventually reach endingNode
  // make the shortest path undefined.
  std::vector<bool> affected(n, false);

  for (size_t u = 0; u < n; ++u) {
    if (!std::isfinite(dist[u]))
      continue;

    for (size_t i = 0; i < g.adjacency[u].size(); ++i) {
      size_t v = g.adjacency[u][i];
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
      size_t v = g.adjacency[u][i];

      if (!affected[v]) {
        affected[v] = true;
        q.push(v);
      }
    }
  }

  // A reachable negative cycle that can reach endingNode
  // means there is no finite shortest path.
  for (size_t u = 0; u < n; ++u) {
    if (affected[u] && canReachEnd[u])
      throw std::runtime_error(
          "No finite shortest path exists: reachable negative cycle detected.");
  }

  // Handle the special case where startingNode == endingNode.
  // Find the cheapest non-empty cycle returning to startingNode.
  if (startingNode == endingNode) {
    Weight bestDistance = INF;
    size_t bestNode = n;

    for (size_t u = 0; u < n; ++u) {
      if (!std::isfinite(dist[u]))
        continue;

      for (size_t i = 0; i < g.adjacency[u].size(); ++i) {
        size_t v = g.adjacency[u][i];

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
      throw std::runtime_error(
          "No cycle returning to the starting node exists.");

    // Reconstruct the cycle
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

  // Check whether endingNode is reachable
  if (!std::isfinite(dist[endingNode]))
    throw std::runtime_error("No path exists between the specified nodes.");

  // Reconstruct the shortest path
  std::vector<size_t> path;
  size_t current = endingNode;

  while (current != startingNode) {
    if (current >= n || parent[current] == n)
      throw std::runtime_error("Failed to reconstruct the shortest path.");

    path.push_back(current);
    current = parent[current];
  }

  path.push_back(startingNode);

  std::reverse(path.begin(), path.end());

  return path;
}

void Graph::exportDot(const std::string &filename) const {
  std::ofstream file(filename);

  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + filename);
  }

  file << "digraph G {\n";
  file << "    node [shape=circle];\n";

  // Export nodes
  for (size_t i = 0; i < adjacency.size(); ++i) {
    file << "    " << i << ";\n";
  }

  // Export edges
  for (size_t from = 0; from < adjacency.size(); ++from) {
    for (size_t to : adjacency[from]) {
      file << "    " << from << " -> " << to << ";\n";
    }
  }

  file << "}\n";

  if (!file.good()) {
    throw std::runtime_error("Failed to write DOT file: " + filename);
  }
}
