#include "shape_functions.hpp"
#include <Eigen/Core>
#include <Eigen/Eigen>

double computeVolume(TetrahedralMesh &mesh, size_t tetId) {
  Tetrahedron &currentTet = mesh.tets[tetId];

  Vertex &p0 = mesh.vertices[currentTet[0]];
  Vertex &p1 = mesh.vertices[currentTet[1]];
  Vertex &p2 = mesh.vertices[currentTet[2]];
  Vertex &p3 = mesh.vertices[currentTet[3]];

  Eigen::Vector3d v1(p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]);
  Eigen::Vector3d v2(p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]);
  Eigen::Vector3d v3(p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]);

  Eigen::Matrix3d mat;
  mat.col(0) = v1;
  mat.col(1) = v2;
  mat.col(2) = v3;

  double volume = std::abs(mat.determinant()) / 6.0;

  return volume;
}

DiffTetShapeFunctions computeDiffShapeFunctions(TetrahedralMesh &mesh,
                                                size_t tetId) {
  const auto &tet = mesh.tets[tetId];
  const auto &p0 = mesh.vertices[tet[0]];

  Eigen::Matrix3d D;
  for (int j = 1; j < 4; ++j) {
    const auto &p = mesh.vertices[tet[j]];
    for (int k = 0; k < 3; ++k) {
      D(k, j - 1) = p[k] - p0[k];
    }
  }

  const double det = D.determinant();
  if (std::abs(det) < 1e-12) {
    throw std::runtime_error("Degenerate tetrahedron");
  }

  const Eigen::Matrix3d invT = D.inverse().transpose();

  DiffTetShapeFunctions res;
  res.g1 = invT.col(0);
  res.g2 = invT.col(1);
  res.g3 = invT.col(2);
  res.g0 = -(res.g1 + res.g2 + res.g3);

  return res;
}

Eigen::Matrix3d computeHexaJacobian(const std::array<Vertex, 8> &nodes,
                                    const Vertex &refPoint) {
  const double xi = refPoint.x();
  const double eta = refPoint.y();
  const double zeta = refPoint.z();

  const int signs[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                           {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};

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
