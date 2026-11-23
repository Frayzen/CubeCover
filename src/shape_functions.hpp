#ifndef SHAPE_FUNCTIONS_HPP
#define SHAPE_FUNCTIONS_HPP

#include "constants.hpp"
#include "ids.hpp"
#include <Eigen/Core>

Eigen::Matrix3d computeHexaJacobian(const std::array<Vertex, 8> &nodes,
                                    const Vertex &refPoint) {
  const double xi = refPoint.x();
  const double eta = refPoint.y();
  const double zeta = refPoint.z();

  const int signs[8][3] = {{-1, -1, -1}, {-1, -1, 1}, {-1, 1, 1}, {-1, 1, -1},
                           {1, -1, -1},  {1, -1, 1},  {1, 1, 1},  {1, 1, -1}};

  Eigen::Matrix3d J = Eigen::Matrix3d::Zero();

  for (int i = 0; i < 8; ++i) {
    const double a = signs[i][0];
    const double b = signs[i][1];
    const double c = signs[i][2];

    const double dN_dxi =
        a * (b > 0 ? eta : 1.0 - eta) * (c > 0 ? zeta : 1.0 - zeta);

    const double dN_deta =
        b * (a > 0 ? xi : 1.0 - xi) * (c > 0 ? zeta : 1.0 - zeta);

    const double dN_dzeta =
        c * (a > 0 ? xi : 1.0 - xi) * (b > 0 ? eta : 1.0 - eta);

    Eigen::Vector3d dN(dN_dxi, dN_deta, dN_dzeta);

    J += nodes[i] * dN.transpose();
  }

  return J;
}

struct DiffTetShapeFunctions {
  Vertex g1;
  Vertex g2;
  Vertex g3;
  Vertex g4;
};

inline DiffTetShapeFunctions
computeDiffShapeFunctions(std::array<Vertex, NB_NODE_PER_HEX> hex) {
  DiffTetShapeFunctions res;

  // Tetrahedron reference coordinates
  const Vertex refPoint(0.25, 0.25, 0.25);

  const Eigen::Matrix3d J = computeHexaJacobian(hex, refPoint);

  const Eigen::Matrix3d invJ = J.inverse();

  // Reference shape function gradients
  const Eigen::Vector3d dN1(-1.0, -1.0, -1.0);
  const Eigen::Vector3d dN2(1.0, 0.0, 0.0);
  const Eigen::Vector3d dN3(0.0, 1.0, 0.0);
  const Eigen::Vector3d dN4(0.0, 0.0, 1.0);

  // Transform gradients to physical coordinates
  res.g1 = invJ.transpose() * dN1;
  res.g2 = invJ.transpose() * dN2;
  res.g3 = invJ.transpose() * dN3;
  res.g4 = invJ.transpose() * dN4;

  return res;
}
#endif
