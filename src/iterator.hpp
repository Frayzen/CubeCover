#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "constants.hpp"
#include "mesh_loader.hpp"
#include <Eigen/SparseCore>
class CubeCoverIterator {

private:
  size_t getRow(size_t tetId, size_t localNode, size_t coord);

  void fillLaplaceTop();
  void fillConstraintTop();

public:
  CubeCoverIterator(TetrahedralMesh &m, std::vector<Eigen::Matrix3d> &frames)
      : mesh(m), frames(frames),
        sizeLaplace(mesh.tets.size() * NB_UNKNOWN_PER_TET),
        sizeConstraint(mesh.amountInternalFace * NB_NODE_PER_FACE * DIMENSION) {
    size_t N = sizeLaplace + sizeConstraint;
    A = Eigen::SparseMatrix<double>(N, N);
    A.reserve(sizeLaplace + sizeConstraint);
    b = Eigen::VectorXd::Zero(N);
    x = Eigen::VectorXd::Zero(N);
  };

  void initialize();

  void solve();

  TetrahedralMesh &mesh;
  std::vector<Eigen::Matrix3d> frames;
  size_t sizeConstraint;
  size_t sizeLaplace;
  Eigen::SparseMatrix<double> A;
  Eigen::VectorX<double> b;
  Eigen::VectorX<double> x;
};

#endif
