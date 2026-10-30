#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>

using AdjacencyWeight = std::vector<std::vector<double>>;

class Graph {
public:
  std::vector<std::vector<std::pair<size_t, size_t>>> adjacency;

  Graph(size_t amountNode) : adjacency(amountNode) {};

  size_t nbNode() { return adjacency.size(); };
  double getWeight(AdjacencyWeight &w, size_t from, size_t to);
};

std::vector<size_t> shortestPath(Graph &g, AdjacencyWeight &weights,
                                 size_t startingNode);

#endif
