#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <string>
#include <vector>

using AdjacencyWeight = std::vector<std::vector<double>>;

using GraphNodeId = size_t;

class Graph {
public:
  std::vector<std::vector<ssize_t>> adjacency;

  Graph(size_t amountNode) : adjacency(amountNode) {};

  size_t nbNode() { return adjacency.size(); };
  double getWeight(AdjacencyWeight &w, size_t from, size_t to);
  GraphNodeId createNode() {
    adjacency.push_back({});
    return adjacency.size() - 1;
  }

  void exportDot(const std::string &filename) const;
};

std::vector<size_t> bellmanFordShortestPath(Graph &g, AdjacencyWeight &weights,
                                            size_t startingNode,
                                            size_t endingNode);

#endif
