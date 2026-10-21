#include "trilaterate.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>

Vertex getNormal(const Vertex plane[3]) {
  const Vertex u = plane[1] - plane[0];
  const Vertex v = plane[2] - plane[0];
  const Vertex cross = u.cross(v);

  const double scale = u.norm() * v.norm();
  if (scale == 0.0 || cross.norm() < 1e-12 * scale) {
    throw std::invalid_argument("Degenerate triangle in tipProject");
  }

  return cross.normalized();
}

Vertex tipProject(const Vertex oldFace[3], const Vertex newFace[3],
                  Vertex oldTip) {
  const Vertex nOld = getNormal(oldFace);
  Vertex nNew = getNormal(newFace);

  const double height = (oldTip - oldFace[0]).dot(nOld);
  const Vertex projected = oldTip - height * nOld;

  const Vertex e0 = oldFace[1] - oldFace[0];
  const Vertex e1 = oldFace[2] - oldFace[0];
  const Vertex q = projected - oldFace[0];

  Eigen::Matrix2d G;
  G << e0.dot(e0), e0.dot(e1), e0.dot(e1), e1.dot(e1);

  Eigen::Vector2d rhs(q.dot(e0), q.dot(e1));

  if (std::abs(G.determinant()) < 1e-24) {
    // throw std::invalid_argument("Degenerate old face in tipProject");
    std::cerr << "Degenerate old face in tipProject";
  }

  const Eigen::Vector2d uv = G.ldlt().solve(rhs);

  const double w0 = 1.0 - uv[0] - uv[1];
  const double w1 = uv[0];
  const double w2 = uv[1];

  const Vertex projectedNew =
      w0 * newFace[0] + w1 * newFace[1] + w2 * newFace[2];

  return projectedNew + height * nNew;
}
