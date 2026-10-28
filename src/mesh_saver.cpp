#include "mesh_saver.hpp"

std::vector<Edge> buildEdges(const std::vector<Tetrahedron> &tets) {
  std::vector<Edge> edges;

  for (const auto &t : tets) {
    for (size_t i = 0; i < NB_NODE_PER_TET; ++i) {
      for (size_t j = i + 1; j < NB_NODE_PER_TET; ++j) {
        Edge e;
        e.v[0] = std::min(t[i], t[j]);
        e.v[1] = std::max(t[i], t[j]);
        edges.push_back(e);
      }
    }
  }
  return edges;
}
