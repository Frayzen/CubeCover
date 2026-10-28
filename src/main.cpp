#include "constants.hpp"
#include "edge_graph.hpp"
#include "iterator.hpp"
#include "mesh_loader.hpp"
#include "shape_functions.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <exception>
#include <iostream>

struct BBox {
  Vertex max;
  Vertex min;
};

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <mesh-file.vtu|mesh-file.vtk>\n";
    return EXIT_FAILURE;
  }

  TetrahedralMesh mesh;
  if (!loadTetrahedralMesh(argv[1], mesh)) {
    return EXIT_FAILURE;
  }

  EdgeGraph g(mesh);
  AdjacencyWeight weightX = computeWeights(g, Vertex(1, 0, 0));
  AdjacencyWeight weightY = computeWeights(g, Vertex(0, 1, 0));
  AdjacencyWeight weightZ = computeWeights(g, Vertex(0, 0, 1));

  // In the order
  // Vertex(0, 0, 0), Vertex(0, 0, 1), Vertex(0, 1, 1), Vertex(0, 1, 0),
  // Vertex(1, 0, 0), Vertex(1, 0, 1), Vertex(1, 1, 1), Vertex(1, 1, 0),

  // Shear cube
  std::array<size_t, NB_CUBE_CORNER> cornerIds = {0, 1, 7, 6, 2, 3, 5, 4};

  // Hole cube
  // std::array<size_t, NB_CUBE_CORNER> cornerIds = {3, 15, 14, 2, 7, 11, 10,
  // 6};
  //
  // std::array<size_t, NB_CUBE_CORNER> = {3, 15, 14, 1, 7, 11, 10, 9};

  std::array<Vertex, NB_CUBE_CORNER> corners;
  for (size_t i = 0; i < NB_CUBE_CORNER; i++)
    corners[i] = mesh.vertices[cornerIds[i]];

  // FRAME FIELD

  // Eigen::Matrix3d frame;
  // Shear cube
  // frame.col(0) << 1, 0, -0.2;
  // frame.col(1) << 0, 1, -0.2;
  // frame.col(2) << 0, 0, 1;

  // Hole cube
  // frame.col(0) << 1, 0, 0;
  // frame.col(1) << 0, 3 * 0.165, 0;
  // frame.col(2) << 0, 0, 1;

  // for (size_t i = 0; i < DIMENSION; i++)
  //   frame.col(i) /=

  std::vector<Eigen::Matrix3d> frames(mesh.tets.size());

  for (size_t t = 0; t < mesh.tets.size(); t++) {
    auto &tet = mesh.tets[t];
    Eigen::Matrix3d jacobian =
        computeHexaJacobian(corners, tet.getBarycenter(mesh.vertices));
    // frames[t] = jacobian.transpose().inverse();
    frames[t] = Eigen::Matrix3d::Identity();
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

  CubeCoverIterator system(mesh, frames, cornerIds);
  system.initialize();

  std::vector<std::string> meshFiles;

  try {
    system.iterate();
  } catch (std::exception &e) {
    std::cout << "Error " << e.what() << std::endl;
    exit(1);
  }

  std::string filename = "output.vtu";
  system.saveMesh(filename);

  return EXIT_SUCCESS;
}
