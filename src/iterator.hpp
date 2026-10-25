#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "constants.hpp"
#include "mesh_loader.hpp"
#include <Eigen/SparseCore>
class CubeCoverIterator {

private:
  size_t getId(size_t tetId, size_t localNode, size_t coord);

  void fillLaplaceTop();
  void debugSolvedGradients();
  void fillConstraintTop();
  void applySpanTree();
  void solve();
  void buildCoordsToRemove();
  void removeOneCoord(size_t coordId);
  Vertex getVertex(size_t tetId, size_t localNodeId);
  std::array<Vertex, NB_NODE_PER_TET> getTetrahedraVertices(size_t tetId);

public:
  CubeCoverIterator(TetrahedralMesh &m, std::vector<Eigen::Matrix3d> &frames)
      : mesh(m), frames(frames), sizeLaplace(mesh.vertices.size() * DIMENSION),
        sizeConstraint(mesh.amountInternalFace * NB_NODE_PER_FACE * DIMENSION) {
    const size_t N = sizeLaplace + sizeConstraint;
    // const size_t N = sizeLaplace;

    A = Eigen::SparseMatrix<double>(N, N);
    A.reserve(N);

    b = Eigen::VectorXd::Zero(N);
    x = Eigen::VectorXd::Zero(N);

    scalarLocation.resize(sizeLaplace);

    for (size_t i = 0; i < sizeLaplace; ++i)
      scalarLocation[i] = static_cast<int>(i);
  }

  void initialize();
  void iterate();
  void saveMesh(std::string path);

  TetrahedralMesh &mesh;
  std::vector<Eigen::Matrix3d> frames;
  size_t sizeConstraint;
  size_t sizeLaplace;
  Eigen::SparseMatrix<double> A;
  Eigen::VectorX<double> b;
  Eigen::VectorX<double> x;

  std::vector<size_t> spanTree;

  Eigen::VectorX<double> spannedX;
  std::vector<double> fixedValues;

  // negative means fixedValues,
  // positive means x (or spannedX)
  std::vector<int> scalarLocation;
  std::vector<int> coordsToRemove;
};

#endif
