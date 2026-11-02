#include "constants.hpp"
#include "edge_graph.hpp"
#include "graph.hpp"
#include "iterator.hpp"
#include "mesh_loader.hpp"
#include "mesh_saver.hpp"
#include "path_saver.hpp"
#include "region_divider.hpp"
#include "shape_functions.hpp"
#include "signed_axis.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <cstdlib>
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

  EdgeGraph edgeGraph(mesh);

  AdjacencyWeight weights[] = {
      computeEdgeWeights(edgeGraph, Vertex(1, 0, 0)),
      computeEdgeWeights(edgeGraph, Vertex(0, 1, 0)),
      computeEdgeWeights(edgeGraph, Vertex(0, 0, 1)),
  };

  RegionDivider regionGraph(edgeGraph);

  double scores[DIMENSION] = {0, 0, 0};
  size_t startEdge[DIMENSION] = {0, 0, 0};
  for (size_t i = 0; i < mesh.nbBoundaryEdge; i++) {
    Edge &e = mesh.edges[i];
    Vertex edgeDir = e.computeDir(mesh);
    for (size_t d = 0; d < DIMENSION; d++) {
      Vertex axis = Vertex::Zero();
      axis[d] = 1;
      double score = abs(axis.dot(edgeDir));
      if (score > scores[d]) {
        scores[d] = score;
        startEdge[d] = i;
      }
    }
  }
  size_t i = 0;
  for (size_t axis = 0; axis < DIMENSION; axis++) {
    size_t startNode = startEdge[axis];
    std::vector<size_t> res =
        bellmanFordShortestPath(edgeGraph, weights[axis], startNode, startNode);
    std::vector<bool> edgeIncluded(mesh.edges.size(), false);
    for (size_t i = 1; i < res.size(); i++) {
      size_t e = res[i];
      size_t next = res[i - 1];
      edgeIncluded[e] = true;
    }
    regionGraph.insertLoop(res, axis);
  }
  std::vector<Vertex> centerEdges(mesh.edges.size());
  for (size_t i = 0; i < mesh.edges.size(); i++)
    centerEdges[i] = mesh.edges[i].computeCenter(mesh);

  savePaths("paths.vtu", regionGraph.paths, centerEdges);

  regionGraph.saveMesh("regions.vtu");

  std::vector<NodeId> corners = regionGraph.computeCorners();
  std::vector<Path> paths = regionGraph.computePaths(corners);
  std::vector<PerDimNormals> colorization = regionGraph.colorizeMesh(paths);

  // In the order
  // Vertex(0, 0, 0), Vertex(0, 0, 1), Vertex(0, 1, 1), Vertex(0, 1, 0),
  // Vertex(1, 0, 0), Vertex(1, 0, 1), Vertex(1, 1, 1), Vertex(1, 1, 0),

  PerDimNormals cornerNormalsOrder[NB_NODE_PER_HEX] = {
      {NEGATIVE, NEGATIVE, NEGATIVE}, {NEGATIVE, NEGATIVE, POSITIVE},
      {NEGATIVE, POSITIVE, POSITIVE}, {NEGATIVE, POSITIVE, NEGATIVE},
      {POSITIVE, NEGATIVE, NEGATIVE}, {POSITIVE, NEGATIVE, POSITIVE},
      {POSITIVE, POSITIVE, POSITIVE}, {POSITIVE, POSITIVE, NEGATIVE},
  };

  // Shear cube
  // std::array<size_t, NB_NODE_PER_HEX> cornerIds = {0, 1, 7, 6, 2, 3, 5, 4};

  // Hole cube
  // std::array<size_t, NB_CUBE_CORNER> cornerIds = {3, 15, 14, 2, 7, 11, 10,
  // 6};
  //
  // std::array<size_t, NB_CUBE_CORNER> = {3, 15, 14, 1, 7, 11, 10, 9};

  // // FRAME FIELD

  std::vector<Eigen::Matrix3d> frames(mesh.tets.size());

  std::array<Vertex, NB_NODE_PER_HEX> cornerVertices;
  std::array<NodeId, NB_NODE_PER_HEX> cornerIds;

  for (size_t i = 0; i < NB_NODE_PER_HEX; i++) {
    PerDimNormals target = cornerNormalsOrder[i];
    for (size_t j = 0; j < NB_NODE_PER_HEX; j++) {
      RegionId regionId = regionGraph.getNodeRegion(corners[j]);
      PerDimNormals normals = regionGraph.getRegionNormals(regionId);
      if (normals == target) {
        cornerVertices[i] = mesh.vertices[corners[j]];
        cornerIds[i] = corners[j];
        break;
      }
    }
  }

  for (size_t t = 0; t < mesh.tets.size(); t++) {
    auto &tet = mesh.tets[t];
    Eigen::Matrix3d jacobian =
        computeHexaJacobian(cornerVertices, tet.getBarycenter(mesh.vertices));
    frames[t] = jacobian.inverse();
    // frames[t] = Eigen::Matrix3d::Identity();
  }
  // for (Face &f : mesh.faces) {
  //   if (f.adjacentCells[1] != UNDEFINED_ELEMENT)
  //     continue; // NOT BOUNDARY
  //   size_t cell = f.adjacentCells[0];
  //   frames[cell] = Eigen::Matrix3d::Identity();
  // }

  // std::cout << "Loaded tetrahedral mesh successfully.\n"
  //           << "Vertices: " << mesh.vertices.size() << '\n'
  //           << "Tetrahedra: " << mesh.tets.size() << '\n';

  // BBox bbox = {
  //     .max = {mind, mind, mind},
  //     .min = {maxd, maxd, maxd},
  // };

  // for (size_t i = 0; i < mesh.vertices.size(); i++) {
  //   auto &v = mesh.vertices[i];
  //   for (size_t coord = 0; coord < DIMENSION; coord++) {
  //     if (v[coord] < bbox.min[coord])
  //       bbox.min[coord] = v[coord];
  //     if (v[coord] > bbox.max[coord])
  //       bbox.max[coord] = v[coord];
  //   }
  // }

  CubeCoverIterator system(mesh, frames, colorization);
  system.initialize();

  std::vector<std::string> meshFiles;

  try {
    system.solve();
  } catch (std::exception &e) {
    std::cout << "Error " << e.what() << std::endl;
    exit(1);
  }
  system.saveMesh("polycube.vtu");

  // std::string filename = "output.vtu";
  // system.saveMesh(filename);

  return EXIT_SUCCESS;
}
