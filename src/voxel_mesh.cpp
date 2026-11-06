#include "voxel_mesh.hpp"
#include "constants.hpp"
#include "signed_axis.hpp"
#include "slice_manager.hpp"
#include "vertex.hpp"
#include <Eigen/Eigen>
#include <algorithm>
#include <cmath>

static bool triangleBoxIntersection(const Vertex &A, const Vertex &B,
                                    const Vertex &C, const Vertex &boxMin,
                                    const Vertex &boxMax) {
  const Vertex center = (boxMin + boxMax) * 0.5;
  const Vertex halfSize = (boxMax - boxMin) * 0.5;

  // Translate triangle to box-centered coordinates.
  const Vertex v0 = A - center;
  const Vertex v1 = B - center;
  const Vertex v2 = C - center;

  const Vertex edges[3] = {v1 - v0, v2 - v1, v0 - v2};

  const Vertex boxAxes[3] = {Vertex::UnitX(), Vertex::UnitY(), Vertex::UnitZ()};

  // Test an axis for separation.
  auto separated = [&](const Vertex &axis) -> bool {
    if (axis.squaredNorm() < 1e-24)
      return false;

    const double p0 = v0.dot(axis);
    const double p1 = v1.dot(axis);
    const double p2 = v2.dot(axis);

    const double triMin = std::min({p0, p1, p2});
    const double triMax = std::max({p0, p1, p2});

    const double radius = halfSize.dot(axis.cwiseAbs());

    return triMin > radius || triMax < -radius;
  };

  // 1. Box coordinate axes.
  for (const Vertex &axis : boxAxes) {
    if (separated(axis))
      return false;
  }

  // 2. Triangle normal.
  const Vertex normal = edges[0].cross(edges[1]);

  if (separated(normal))
    return false;

  // 3. Cross products of triangle edges and box axes.
  for (const Vertex &edge : edges) {
    for (const Vertex &axis : boxAxes) {
      if (separated(edge.cross(axis)))
        return false;
    }
  }

  return true;
}

static std::vector<Voxel> triangleVoxels(const Vertex &A, const Vertex &B,
                                         const Vertex &C, double voxelSize,
                                         const Vertex &origin) {
  if (voxelSize <= 0.0)
    throw std::invalid_argument("voxelSize must be positive");

  // Triangle bounding box.
  const Vertex triMin = A.cwiseMin(B).cwiseMin(C);

  const Vertex triMax = A.cwiseMax(B).cwiseMax(C);

  // Convert to integer grid coordinates.
  const Voxel gridMin =
      ((triMin - origin) / voxelSize).array().floor().cast<int>();

  const Voxel gridMax =
      ((triMax - origin) / voxelSize).array().floor().cast<int>();

  std::vector<Voxel> result;

  // Enumerate candidate voxels.
  for (int i = gridMin.x(); i <= gridMax.x(); ++i) {
    for (int j = gridMin.y(); j <= gridMax.y(); ++j) {
      for (int k = gridMin.z(); k <= gridMax.z(); ++k) {
        const Vertex boxMin = origin + voxelSize * Vertex(i, j, k);
        const Vertex boxMax = boxMin + Vertex::Constant(voxelSize);
        if (triangleBoxIntersection(A, B, C, boxMin, boxMax))
          result.emplace_back(i, j, k);
      }
    }
  }

  return result;
}

static void updateNormals(Axis axis, Region &oldRegion, Region &newRegion) {
  IntegerSign oldNormal = oldRegion.signedNormals[axis];
  if (oldNormal == UNDEFINED) {
    oldRegion.signedNormals[axis] = NEGATIVE;
    newRegion.signedNormals[axis] = POSITIVE;
  } else if (oldNormal == ZERO)
    newRegion.signedNormals[axis] = ZERO;
  else if (oldNormal == NEGATIVE)
    newRegion.signedNormals[axis] = ZERO;
  else
    oldRegion.signedNormals[axis] = ZERO;
}

void VoxelMesh::slice(Axis axis, SliceId splittingSlice) {
  auto &sliceManager = sliceManagers[axis];
  std::vector<Region> createdRegions;
  for (Region &region : regions) {
    assert(region.startSlices[axis] != splittingSlice);
    if (region.containsSlice(axis, splittingSlice)) {
      Region newRegion(region);
      newRegion.startSlices[axis] = splittingSlice;
      region.endSlices[axis] = splittingSlice;
      std::vector<VoxelId> newRegionIds;
      std::vector<VoxelId> oldRegionIds;
      for (VoxelId voxelId : region) { // new region is above
        if (voxels[voxelId][axis] > splittingSlice)
          newRegionIds.push_back(voxelId);
        else
          oldRegionIds.push_back(voxelId);
      }
      region.voxels = oldRegionIds;
      newRegion.voxels = newRegionIds;
      updateNormals(axis, region, newRegion);
      createdRegions.push_back(newRegion);
    }
  }
  regions.insert(regions.begin(), createdRegions.begin(), createdRegions.end());
}

static void deduplicateVoxels(std::vector<Voxel> &voxels) {
  std::sort(voxels.begin(), voxels.end(), [](const Voxel &a, const Voxel &b) {
    if (a.x() != b.x())
      return a.x() < b.x();
    if (a.y() != b.y())
      return a.y() < b.y();
    return a.z() < b.z();
  });
  voxels.erase(std::unique(voxels.begin(), voxels.end()), voxels.end());
}

void VoxelMesh::buildMesh() {

  for (const Vertex &v : tetMesh.vertices)
    bbox.update(v);

  origin = bbox.min;

  const Vertex diff = bbox.max - bbox.min;

  voxelSize = diff.maxCoeff() / 50.0;

  if (voxelSize <= 0.0)
    throw std::runtime_error("Invalid voxel size.");

  for (size_t faceId = 0; faceId < tetMesh.nbBoundaryFace; faceId++) {
    const Face &face = tetMesh.faces[faceId];
    const Vertex &v1 = tetMesh.vertices[face[0]];
    const Vertex &v2 = tetMesh.vertices[face[1]];
    const Vertex &v3 = tetMesh.vertices[face[2]];
    auto found = triangleVoxels(v1, v2, v3, voxelSize, origin);

    voxels.insert(voxels.end(), found.begin(), found.end());
  }
  deduplicateVoxels(voxels);
  for (SliceManager &map : sliceManagers)
    map.build(*this);

  Region initialRegion(*this);
  initialRegion.signedNormals = {UNDEFINED, UNDEFINED, UNDEFINED};
  for (Axis d = 0; d < DIMENSION; d++) {
    initialRegion.startSlices[d] = 0;
    initialRegion.endSlices[d] = sliceManagers[d].nbSlice();
  }
  for (VoxelId id = 0; id < voxels.size(); id++)
    initialRegion.voxels.push_back(id);
  regions.push_back(initialRegion);
}
