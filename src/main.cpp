#include "constants.hpp"
#include "edge_graph.hpp"
#include "graph.hpp"
#include "mesh_loader.hpp"
#include "path_saver.hpp"
#include "region_divider.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <cstdlib>
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

  RegionDivider regionDivider(edgeGraph);

  double scores[DIMENSION] = {0, 0, 0};
  EdgeId startEdge[DIMENSION] = {0, 0, 0};
  for (EdgeId edgeId = 0; edgeId < mesh.nbBoundaryEdge; edgeId++) {
    Edge &e = mesh.edges[edgeId];
    Vertex edgeDir = e.computeDir(mesh);
    for (size_t d = 0; d < DIMENSION; d++) {
      Vertex axis = Vertex::Zero();
      axis[d] = 1;
      double score = abs(axis.dot(edgeDir));
      if (score > scores[d]) {
        scores[d] = score;
        startEdge[d] = edgeId;
      }
    }
  }
  auto addLoops = [&]() {
    for (size_t axis = 0; axis < DIMENSION; axis++) {
      size_t startNode = startEdge[axis];
      std::vector<size_t> res = bellmanFordShortestPath(
          edgeGraph, weights[axis], startNode, startNode);
      regionDivider.insertLoop(res, axis);
    }
  };
  addLoops();
  addLoops();
  std::vector<Vertex> centerEdges(mesh.edges.size());
  for (size_t i = 0; i < mesh.edges.size(); i++)
    centerEdges[i] = mesh.edges[i].computeCenter(mesh);

  double quality = regionDivider.measureQuality();

  std::srand(1234);
  RegionDivider tmp = regionDivider;
  for (size_t i = 0; i < 1; i++) {
    NodeId startNode = std::rand() % mesh.nbBoundaryNode;
    Axis randomAxis = std::rand() % DIMENSION;
    std::vector<size_t> res = bellmanFordShortestPath(
        edgeGraph, weights[randomAxis], startNode, startNode);
    tmp.insertLoop(res, randomAxis);
  }
  savePaths("paths.vtu", tmp.paths, centerEdges);

  // savePaths("paths.vtu", regionDivider.paths, centerEdges);

  //   regionGraph.saveMesh("regions.vtu");

  //   std::vector<NodeId> corners = regionGraph.computeCorners();
  //   std::vector<Path> paths = regionGraph.computePaths(corners);
  //   std::vector<PerDimNormals> colorization =
  //   regionGraph.colorizeMesh(paths);

  //   PerDimNormals cornerNormalsOrder[NB_NODE_PER_HEX] = {
  //       {NEGATIVE, NEGATIVE, NEGATIVE}, {NEGATIVE, NEGATIVE, POSITIVE},
  //       {NEGATIVE, POSITIVE, POSITIVE}, {NEGATIVE, POSITIVE, NEGATIVE},
  //       {POSITIVE, NEGATIVE, NEGATIVE}, {POSITIVE, NEGATIVE, POSITIVE},
  //       {POSITIVE, POSITIVE, POSITIVE}, {POSITIVE, POSITIVE, NEGATIVE},
  //   };

  //   // FRAME FIELD

  //   std::vector<Eigen::Matrix3d> frames(mesh.tets.size());

  //   std::array<Vertex, NB_NODE_PER_HEX> cornerVertices;
  //   std::array<NodeId, NB_NODE_PER_HEX> cornerIds;

  //   for (size_t i = 0; i < NB_NODE_PER_HEX; i++) {
  //     PerDimNormals target = cornerNormalsOrder[i];
  //     for (size_t j = 0; j < NB_NODE_PER_HEX; j++) {
  //       RegionId regionId = regionGraph.getNodeRegion(corners[j]);
  //       PerDimNormals normals = regionGraph.getRegionNormals(regionId);
  //       if (normals == target) {
  //         cornerVertices[i] = mesh.vertices[corners[j]];
  //         cornerIds[i] = corners[j];
  //         break;
  //       }
  //     }
  //   }

  //   for (size_t t = 0; t < mesh.tets.size(); t++) {
  //     auto &tet = mesh.tets[t];
  //     Eigen::Matrix3d jacobian =
  //         computeHexaJacobian(cornerVertices,
  //         tet.getBarycenter(mesh.vertices));
  //     frames[t] = jacobian.inverse();
  //   }

  //   CubeCoverIterator system(mesh, frames, colorization);
  //   system.initialize();

  //   std::vector<std::string> meshFiles;

  //   try {
  //     system.solve();
  //   } catch (std::exception &e) {
  //     std::cout << "Error " << e.what() << std::endl;
  //     exit(1);
  //   }
  //   system.saveMesh("polycube.vtu");

  //   return EXIT_SUCCESS;
}
