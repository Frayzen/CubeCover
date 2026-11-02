#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "constants.hpp"
#include "mesh_loader.hpp"
#include "region_divider.hpp"
#include <Eigen/SparseCore>
class CubeCoverIterator {

private:
  size_t getId(size_t tetId, size_t localNode, size_t coord);

  void fillLaplaceTop();
  void debugSolvedGradients();
  void fillConstraints();
  Vertex getVertex(size_t tetId, size_t localNodeId);
  std::array<Vertex, NB_NODE_PER_TET> getTetrahedraVertices(size_t tetId);

public:
  CubeCoverIterator(TetrahedralMesh &m, std::vector<Eigen::Matrix3d> &frames,
                    std::vector<PerDimNormals> nodeToAxisConstraints)
      : mesh(m), frames(frames), sizeLaplace(mesh.vertices.size() * DIMENSION),
        sizeConstraint(0), nodeToAxisConstraints(nodeToAxisConstraints) {

    for (NodeId nodeId = 0; nodeId < mesh.nbBoundaryNode; nodeId++) {
      for (Axis d = 0; d < DIMENSION; d++) {
        assert(nodeToAxisConstraints[nodeId][d] != UNDEFINED);
        if (nodeToAxisConstraints[nodeId][d] != ZERO)
          sizeConstraint++;
      }
    }

    const size_t N = sizeLaplace + sizeConstraint;

    A = Eigen::SparseMatrix<double>(N, N);
    A.reserve(N);

    b = Eigen::VectorXd::Zero(N);
    x = Eigen::VectorXd::Zero(N);
    for (size_t i = 0; i < mesh.vertices.size(); i++) {
      for (size_t d = 0; d < DIMENSION; d++)
        x[DIMENSION * i + d] = mesh.vertices[i][d];
    }
  }

  void initialize();
  void solve();
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

  std::vector<PerDimNormals> nodeToAxisConstraints;
};

#endif
