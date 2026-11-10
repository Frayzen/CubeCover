#ifndef JUNCTIONS_HPP
#define JUNCTIONS_HPP
#include "ids.hpp"
#include "mesh_loader.hpp"

struct SurfacePatches {
  std::vector<int> faceToPatch;
  std::vector<std::vector<FaceId>> patches;
};

SurfacePatches buildSurfacePatches(TetrahedralMesh &mesh);
std::vector<TetId> computeMedialAxisJunctions(const TetrahedralMesh &tetMesh,
                                              const SurfacePatches &surface);

std::vector<Vertex>
computeJunctionClusters(TetrahedralMesh &mesh,
                        const std::vector<TetId> &junctionTetIds);

#endif
