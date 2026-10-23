#ifndef SHAPE_FUNCTIONS_HPP
#define SHAPE_FUNCTIONS_HPP

#include "mesh_loader.hpp"
#include <iostream>
struct DiffTetShapeFunctions {
  Vertex g0;
  Vertex g1;
  Vertex g2;
  Vertex g3;

  inline Vertex &operator[](size_t i) {
    switch (i) {
    case 0:
      return g0;
    case 1:
      return g1;
    case 2:
      return g2;
    case 3:
      return g3;
    default:
      std::cerr << "Wrong value provided " << i << std::endl;
      assert(false);
    }
  }

  inline const Vertex &operator[](size_t i) const {
    switch (i) {
    case 0:
      return g0;
    case 1:
      return g1;
    case 2:
      return g2;
    case 3:
      return g3;
    default:
      std::cerr << "Wrong value provided " << i << std::endl;
      assert(false);
    }
  }
};

double computeVolume(TetrahedralMesh &mesh, size_t tetId);

DiffTetShapeFunctions computeDiffShapeFunctions(TetrahedralMesh &mesh,
                                                size_t tetId);
#endif
