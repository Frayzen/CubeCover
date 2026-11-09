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

  VoxelMesh mesh(tetMesh);
  for (Axis d = 0; d < DIMENSION; d++) {
    size_t nbSlice = mesh.sliceManagers[d].nbSlice();
    mesh.cut(d, nbSlice / 2);
  }

  for (Axis d = 0; d < DIMENSION; d++)
    mesh.regionGraphs[d].exportDot("graph_dim" + std::to_string(d) + ".dot");

  // mesh.splitByConnectedRegions();
  std::vector<VoxelId> corners = mesh.findCorners();
  std::vector<std::vector<VoxelId>> paths = mesh.findPaths();

  VoxelMeshSaver saver(mesh);

  std::vector<int> pathIds(mesh.voxels.size());
  std::fill(pathIds.begin(), pathIds.end(), 0);
  int pathId = 1;
  for (auto path : paths) {
    for (VoxelId v : path)
      pathIds[v] = pathId;
    pathId++;
  }

  std::vector<bool> cornersFlag(mesh.voxels.size());
  std::fill(cornersFlag.begin(), cornersFlag.end(), false);

  std::vector<size_t> colors(mesh.voxels.size());
  std::vector<double> xNormal(mesh.voxels.size());
  std::vector<double> yNormal(mesh.voxels.size());
  std::vector<double> zNormal(mesh.voxels.size());

  std::vector<double> xLayer(mesh.voxels.size());
  std::vector<double> yLayer(mesh.voxels.size());
  std::vector<double> zLayer(mesh.voxels.size());
  size_t regionId = 0;
  for (Region region : mesh.regions) {
    for (VoxelId id : region) {
      colors[id] = regionId;
      xNormal[id] = getIntegerSignValue(region.signedNormals[0]);
      yNormal[id] = getIntegerSignValue(region.signedNormals[1]);
      zNormal[id] = getIntegerSignValue(region.signedNormals[2]);

      xLayer[id] = region.layerIds[0];
      yLayer[id] = region.layerIds[1];
      zLayer[id] = region.layerIds[2];
      if (id == corners[regionId])
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
