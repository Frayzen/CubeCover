#ifndef PRIMAL_REGION_HPP
#define PRIMAL_REGION_HPP
#include "constants.hpp"
#include "ids.hpp"
#include "signed_axis.hpp"

class VoxelMesh;
struct PrimalRegion {
  std::array<SurfaceId, NB_FACE_PER_HEX> surfaceIds;
  std::vector<Voxel> voxels;
  /* VTK-LIKE ORDERING
// 0	-	-	-
// 1	+	-	-
// 2	+	+	-
// 3	-	+	-
// 4	-	-	+
// 5	+	-	+
// 6	+	+	+
// 7	-	+	+
   */
  std::array<Voxel, NB_NODE_PER_HEX> corners;
  std::array<DualRegionId, NB_NODE_PER_HEX> dualRegions;

  PrimalRegion() {
    for (size_t i = 0; i < NB_FACE_PER_HEX; i++)
      surfaceIds[i] = NO_SURFACE;
  }
  inline SurfaceId &getSurfaceId(Axis a, IntegerSign sign) {
    return surfaceIds[a * 2 + (sign == POSITIVE ? 1 : 0)];
  }
  Voxel getRegionNormals(DualRegionId dualRegionId) const;
  std::array<DualRegionId, NB_NODE_PER_QUAD>
  getRegionFromFace(Axis a, IntegerSign sign) const;
  void reorderDualRegionsToVTK(const VoxelMesh &mesh);
};
#endif
