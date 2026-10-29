#ifndef AXIS_GRAPH_HPP
#define AXIS_GRAPH_HPP

#include "mesh_loader.hpp"
class AxisGraph {
private:
  TetrahedralMesh &m;

public:
  AxisGraph(TetrahedralMesh &m) : m(m) {};

  void insertLoop() {}
};

#endif
