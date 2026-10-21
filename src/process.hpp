#ifndef CUBE_COVER_ITERATOR_HPP
#define CUBE_COVER_ITERATOR_HPP
#include "mesh_loader.hpp"
#include "span_tree.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>

struct CubeCoverParameters {
  Eigen::SparseMatrix<double> system;
  Eigen::VectorX<double> knowns;
  Eigen::VectorX<double> unspannedSolution;
  Eigen::VectorX<double> result;
  TetrahedralMesh &mesh;
  std::vector<size_t> spanTree;
  std::vector<Eigen::Matrix3d> frames;
  std::vector<int> valueMapper;

  CubeCoverParameters(TetrahedralMesh &mesh, size_t systemSize)
      : mesh(mesh), frames(mesh.tets.size()), system(systemSize, systemSize),
        knowns(systemSize), result(Eigen::VectorX<double>::Zero(systemSize)),
        spanTree(convertToSpanTree(mesh)) {}
};

void cubeCoverIteration(CubeCoverParameters &params, size_t boundaryFaceId);
#endif
