#include "mesh_loader.hpp"
#include "signed_axis.hpp"
#include "voxel_mesh.hpp"
#include "voxel_mesh_saver.hpp"
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <mesh-file.vtu|mesh-file.vtk>\n";
    return EXIT_FAILURE;
  }

  TetrahedralMesh tetMesh;
  if (!loadTetrahedralMesh(argv[1], tetMesh)) {
    return EXIT_FAILURE;
  }

  // SurfacePatches patches = buildSurfacePatches(tetMesh);
  // std::vector<size_t> patchIds(tetMesh.nbBoundaryFace);
  // std::fill(patchIds.begin(), patchIds.end(), 0);
  // size_t curPatch = 0;
  // for (auto patch : patches.patches) {
  //   for (FaceId f : patch) {
  //     patchIds[f] = curPatch;
  //     assert(f < tetMesh.nbBoundaryFace);
  //   }
  //   curPatch++;
  // }

  // std::vector<TetId> junctions = computeMedialAxisJunctions(tetMesh,
  // patches); std::cout << "FOUND " << junctions.size() << " JUNCITONS " <<
  // std::endl; std::vector<bool> junctionFlag(tetMesh.tets.size());
  // std::fill(junctionFlag.begin(), junctionFlag.end(), false);
  // for (TetId tId : junctions)
  //   junctionFlag[tId] = true;

  // MeshSaver tetMeshSaver(tetMesh, "tetMesh.vtu");
  // tetMeshSaver.addField<ElementType::CellType>("isJunction", junctionFlag);
  // tetMeshSaver.addField<ElementType::FaceBoundaryType>("patchIds", patchIds);
  // tetMeshSaver.saveMesh();
  // tetMeshSaver.saveFaces();

  // std::vector<Vertex> clusters = computeJunctionClusters(tetMesh, junctions);
  // saveNodes(clusters, "clusters.vtu");
  // std::cout << "FOUND " << clusters.size() << " CLUSTERS " << std::endl;

  VoxelMesh mesh(tetMesh);

  // for (size_t i = 0; i < clusters.size(); i++) {
  //   for (size_t j = i; j < clusters.size(); j++) {
  //     for (Axis d = 0; d < DIMENSION; d++) {
  //       if (std::abs(clusters[i][d] - clusters[j][d]) < (3 * mesh.voxelSize))
  //         clusters[i][d] = clusters[j][d];
  //     }
  //   }
  // }

  // for (Vertex cutVertex : clusters) {
  //   for (Axis d = 0; d < DIMENSION; d++) {
  //     double coord = cutVertex[d];
  //     double relative = coord - mesh.origin[d];
  //     SliceId sliceId = std::round(relative / mesh.voxelSize);
  //     mesh.cut(d, sliceId);
  //   }
  // }

  for (Axis d = 0; d < DIMENSION; d++) {
    double nbSlice = mesh.voxels.dimension(d);
    mesh.cut(d, nbSlice / 2);
  }

  for (Axis d = 0; d < DIMENSION; d++)
    mesh.regionGraphs[d].exportDot("graph_dim" + std::to_string(d) + ".dot");
  std::vector<int> pathIds(mesh.nbVoxel);

  std::vector<Voxel> corners = mesh.findCorners();
  std::vector<std::vector<Voxel>> paths = mesh.findPaths();

  std::fill(pathIds.begin(), pathIds.end(), 0);
  int pathId = 1;
  for (auto path : paths) {
    for (Voxel v : path)
      pathIds[mesh.getId(v)] = pathId;
    pathId++;
  }

  VoxelMeshSaver saver(mesh);

  std::vector<bool> cornersFlag(mesh.nbVoxel);
  std::fill(cornersFlag.begin(), cornersFlag.end(), false);
  std::vector<size_t> colors(mesh.nbVoxel);
  std::fill(colors.begin(), colors.end(), NO_REGION);

  std::vector<double> xNormal(mesh.nbVoxel);
  std::vector<double> yNormal(mesh.nbVoxel);
  std::vector<double> zNormal(mesh.nbVoxel);

  std::vector<double> xLayer(mesh.nbVoxel);
  std::vector<double> yLayer(mesh.nbVoxel);
  std::vector<double> zLayer(mesh.nbVoxel);
  RegionId regionId = 0;
  for (Region region : mesh.regions) {
    for (Voxel v : region) {
      VoxelId id = mesh.getId(v);
      colors[id] = regionId;
      xNormal[id] = getIntegerSignValue(region.signedNormals[0]);
      yNormal[id] = getIntegerSignValue(region.signedNormals[1]);
      zNormal[id] = getIntegerSignValue(region.signedNormals[2]);

      xLayer[id] = region.layerIds[0];
      yLayer[id] = region.layerIds[1];
      zLayer[id] = region.layerIds[2];
      if (v == corners[regionId])
        cornersFlag[id] = true;
    }
    regionId++;
  }
  saver.addField("region", colors);
  saver.addField("xNormal", xNormal);
  saver.addField("yNormal", yNormal);
  saver.addField("zNormal", zNormal);
  saver.addField("xLayer", xLayer);
  saver.addField("yLayer", yLayer);
  saver.addField("zLayer", zLayer);
  saver.addField("corner", cornersFlag);
  saver.addField("pathIds", pathIds);

  saver.saveMesh();
}
