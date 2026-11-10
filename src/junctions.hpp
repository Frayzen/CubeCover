#ifndef JUNCTIONS_HPP
#define JUNCTIONS_HPP
#include "ids.hpp"
#include "mesh_loader.hpp"
#include "voxel_mesh.hpp"

struct SurfacePatches {
  const TetrahedralMesh &mesh;
  std::vector<int> faceToPatch;
  std::vector<std::vector<FaceId>> patches;
};

SurfacePatches buildSurfacePatches(TetrahedralMesh &mesh);
std::vector<Voxel> computeMedialAxisJunctions(const VoxelMesh &mesh,
                                              const SurfacePatches &surface);

std::vector<Vertex>
computeJunctionClusters(TetrahedralMesh &mesh,
                        const std::vector<TetId> &junctionTetIds);

#endif
