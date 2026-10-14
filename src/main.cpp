#include "constants.hpp"
#include "mesh_loader.hpp"
#include <Eigen/Eigen>
#include <iostream>

struct BBox {
  Vertex max;
  Vertex min;
};

struct DiffTetShapeFunctions {
  Vertex g0;
  Vertex g1;
  Vertex g2;
  Vertex g3;
};

static double computeArea(TetrahedralMesh &mesh, size_t tetId) {
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

static auto computeDiffShapeFunctions(TetrahedralMesh &mesh, size_t tetId) {
  double data[3][3] = {-999};
  Tetrahedron &currentTet = mesh.tets[tetId];
  Vertex &p0 = mesh.vertices[currentTet[0]];
  for (size_t i = 1; i < TET_NODE_COUNT; i++) {
    Vertex &pi = mesh.vertices[currentTet[i]];
    for (size_t coord = 0; coord < DIMENSION; coord++)
      data[i - 1][coord] = pi[coord] - p0[coord];
  }
  Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor>> mat(data[0]);
  auto inv = mat.inverse();

  DiffTetShapeFunctions res{
      .g0 = {},
      .g1 = inv.row(0),
      .g2 = inv.row(1),
      .g3 = inv.row(2),
  };
  res.g0 = Vertex(1, 1, 1) - res.g1 - res.g2 - res.g3;

  // std::cout << "=============================" << std::endl;
  // std::cout << "GOT AREA " << computeArea(mesh, tetId) << std::endl;
  // std::cout << res.g0 << std::endl;
  // std::cout << res.g1 << std::endl;
  // std::cout << res.g2 << std::endl;
  // std::cout << res.g3 << std::endl;
  return res;
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <mesh-file.vtu|mesh-file.vtk>\n";
    return EXIT_FAILURE;
  }

  TetrahedralMesh mesh;
  if (!loadTetrahedralMesh(argv[1], mesh)) {
    return EXIT_FAILURE;
  }

  std::cout << "Loaded tetrahedral mesh successfully.\n"
            << "Vertices: " << mesh.vertices.size() << '\n'
            << "Tetrahedra: " << mesh.tets.size() << '\n';

  BBox bbox = {
      .max = {mind, mind, mind},
      .min = {maxd, maxd, maxd},
  };

  for (size_t i = 0; i < mesh.vertices.size(); i++) {
    auto &v = mesh.vertices[i];
    for (size_t coord = 0; coord < DIMENSION; coord++) {
      if (v[coord] < bbox.min[coord])
        bbox.min[coord] = v[coord];
      if (v[coord] > bbox.max[coord])
        bbox.max[coord] = v[coord];
    }
  }

  for (size_t i = 0; i < mesh.tets.size(); i++) {
    auto shapeFunctions = computeDiffShapeFunctions(mesh, i);
  }

  return EXIT_SUCCESS;
}
