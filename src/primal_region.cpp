#include "primal_region.hpp"
#include "voxel_mesh.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

void PrimalRegion::reorderDualRegionsToVTK(const VoxelMesh &mesh) {
  constexpr size_t N = NB_NODE_PER_HEX;

  constexpr std::array<std::array<int, DIMENSION>, N> vtkOrder = {{
      {{-1, -1, -1}}, // 0
      {{1, -1, -1}},  // 1
      {{1, 1, -1}},   // 2
      {{-1, 1, -1}},  // 3
      {{-1, -1, 1}},  // 4
      {{1, -1, 1}},   // 5
      {{1, 1, 1}},    // 6
      {{-1, 1, 1}}    // 7
  }};

  const auto oldDualRegions = dualRegions;
  const auto oldCorners = corners;

  std::array<std::array<int, DIMENSION>, N> signs{};
  std::array<DualRegionId, N> reorderedDualRegions{};
  std::array<Voxel, N> reorderedCorners{};

  for (size_t i = 0; i < N; ++i) {
    const auto &region = mesh.regions[oldDualRegions[i]];

    for (size_t axis = 0; axis < DIMENSION; ++axis) {
      const auto &children = region.children[axis];

      const bool isPositive = std::any_of(
          children.begin(), children.end(),
          [&oldDualRegions](const DualRegionInterface &interface) {
            return std::find(oldDualRegions.begin(), oldDualRegions.end(),
                             interface.regionId) != oldDualRegions.end();
          });

      signs[i][axis] = isPositive ? 1 : -1;
    }
  }

  // Verify four positive and four negative regions per axis.
  for (size_t axis = 0; axis < DIMENSION; ++axis) {
    size_t positiveCount = 0;

    for (size_t i = 0; i < N; ++i)
      positiveCount += (signs[i][axis] == 1);

    assert(positiveCount == 4);
  }

  // Reorder both arrays according to VTK node numbering.
  std::array<bool, N> assigned{};

  for (size_t i = 0; i < N; ++i) {
    bool found = false;

    for (size_t vtkId = 0; vtkId < N; ++vtkId) {
      if (signs[i] == vtkOrder[vtkId]) {
        assert(!assigned[vtkId]);

        reorderedDualRegions[vtkId] = oldDualRegions[i];
        reorderedCorners[vtkId] = oldCorners[i];

        assigned[vtkId] = true;
        found = true;
        break;
      }
    }

    assert(found);
  }

  for (bool value : assigned)
    assert(value);

  dualRegions = reorderedDualRegions;
  corners = reorderedCorners;
}

std::array<DualRegionId, NB_NODE_PER_QUAD>
PrimalRegion::getRegionFromFace(Axis a, IntegerSign sign) const {
  const std::array<size_t, NB_NODE_PER_QUAD>
      signedAxisToLocalIds[DIMENSION * 2] = {
          // Axis 0 (X) Sign NEGATIVE
          {0, 3, 7, 4},
          // Axis 0 (X) Sign POSITIVE
          {1, 2, 6, 5},

          // Axis 1 (Y) Sign NEGATIVE
          {0, 1, 5, 4},
          // Axis 1 (Y) Sign POSITIVE
          {3, 2, 6, 7},

          // Axis 2 (Z) Sign NEGATIVE
          {0, 1, 2, 3},
          // Axis 2 (Z) Sign POSITIVE
          {4, 5, 6, 7},
      };
  std::array<DualRegionId, NB_NODE_PER_QUAD> res;
  auto localIds = signedAxisToLocalIds[a * 2 + (sign == POSITIVE ? 1 : 0)];
  for (size_t i = 0; i < NB_NODE_PER_QUAD; i++)
    res[i] = dualRegions[localIds[i]];
  return res;
}

Voxel PrimalRegion::getRegionNormals(DualRegionId dualRegionId) const {
  const std::array<Voxel, NB_NODE_PER_HEX> normals = {
      Voxel(-1, -1, -1), Voxel(+1, -1, -1), Voxel(+1, +1, -1),
      Voxel(-1, +1, -1), Voxel(-1, -1, +1), Voxel(+1, -1, +1),
      Voxel(+1, +1, +1), Voxel(-1, +1, +1),
  };
  size_t i = 0;
  for (DualRegionId rId : dualRegions) {
    if (rId == dualRegionId)
      return normals[i];
    i++;
  }
  throw std::runtime_error("Dual region id not in the regions");
};
