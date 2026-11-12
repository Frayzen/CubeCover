#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "ids.hpp"
#include <string>
#include <vector>

using AdjacencyWeight = std::vector<std::vector<double>>;

class Graph {
public:
  std::vector<std::vector<GraphNodeId>> adjacency;

  Graph(size_t amountNode) : adjacency(amountNode) {};

  size_t nbNode() const { return adjacency.size(); };
  double getWeight(AdjacencyWeight &w, size_t from, size_t to);
  GraphNodeId createNode() {
    adjacency.push_back({});
    return adjacency.size() - 1;
  }

  void exportDot(const std::string &filename) const;
  std::vector<std::vector<GraphNodeId>> parents() const;
  const std::vector<GraphNodeId> &operator[](GraphNodeId nodeId) const {
    return adjacency[nodeId];
  }

  std::vector<GraphNodeId> &operator[](GraphNodeId nodeId) {
    return adjacency[nodeId];
  }
};

std::vector<size_t> bellmanFordShortestPath(Graph &g, AdjacencyWeight &weights,
                                            size_t startingNode,
                                            size_t endingNode);

#endif
