#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "constants.hpp"
#include "mesh_loader.hpp"
#include <Eigen/SparseCore>
class CubeCoverIterator {

private:
  size_t getId(size_t tetId, size_t localNode, size_t coord);

  void fillLaplaceTop();
  void fillConstraintTop();
  Vertex getVertex(size_t tetId, size_t localNodeId);
  std::array<Vertex, NB_NODE_PER_TET> getTetrahedraVertices(size_t tetId);

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
    scalarLocation.resize(N);
    for (size_t i = 0; i < N; i++)
      scalarLocation[i] = i;
  };

  void initialize();
  void solve();
  void applySpanTree(std::vector<size_t> spanTree);
  void saveMesh(std::string path);

  TetrahedralMesh &mesh;
  std::vector<Eigen::Matrix3d> frames;
  size_t sizeConstraint;
  size_t sizeLaplace;
  Eigen::SparseMatrix<double> A;
  Eigen::VectorX<double> b;
  Eigen::VectorX<double> x;

  Eigen::VectorX<double> spannedX;
  std::vector<double> fixedValues;

  // negative means fixedValues,
  // positive means x (or spannedX)
  std::vector<int> scalarLocation;
};

#endif
