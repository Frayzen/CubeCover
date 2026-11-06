#include "voxel_mesh.hpp"
#include "constants.hpp"
#include "signed_axis.hpp"
#include "slice_manager.hpp"
#include "vertex.hpp"
#include <Eigen/Eigen>
#include <algorithm>
#include <cmath>

#include <queue>
#include <unsupported/Eigen/CXX11/Tensor>

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
  std::vector<Region> updatedRegions;
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
      updatedRegions.push_back(region);
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

std::vector<Region> VoxelMesh::getConnectedRegions(RegionId rid) {
  std::vector<Region> result;

  const Region &region = regions[rid];

  // Local tensor dimensions.
  const size_t nx = region.endSlices.x() - region.startSlices.x() + 1;
  const size_t ny = region.endSlices.y() - region.startSlices.y() + 1;
  const size_t nz = region.endSlices.z() - region.startSlices.z() + 1;

  // Tensor representing the voxels belonging to this region.
  Eigen::Tensor<bool, 3> tensor(nx, ny, nz);
  tensor.setConstant(false);

  // Store the original voxel IDs.
  Eigen::Tensor<VoxelId, 3> voxelIds(nx, ny, nz);

  // Populate the tensor.
  for (VoxelId voxelId : region.voxels) {
    const auto &voxel = voxels[voxelId];

    const size_t x = voxel.x() - region.startSlices.x();
    const size_t y = voxel.y() - region.startSlices.y();
    const size_t z = voxel.z() - region.startSlices.z();

    tensor(x, y, z) = true;
    voxelIds(x, y, z) = voxelId;
  }

  Eigen::Tensor<bool, 3> visited(nx, ny, nz);
  visited.setConstant(false);

  // 6-connectivity.
  const std::array<Eigen::Vector3<int>, 6> neighbors = {
      Eigen::Vector3<int>(1, 0, 0), Eigen::Vector3<int>(-1, 0, 0),
      Eigen::Vector3<int>(0, 1, 0), Eigen::Vector3<int>(0, -1, 0),
      Eigen::Vector3<int>(0, 0, 1), Eigen::Vector3<int>(0, 0, -1)};

  using Coord = Eigen::Vector3<size_t>;

  // Find connected components.
  for (size_t x = 0; x < nx; ++x) {
    for (size_t y = 0; y < ny; ++y) {
      for (size_t z = 0; z < nz; ++z) {
        if (!tensor(x, y, z) || visited(x, y, z))
          continue;

        // Create a new region for this component.
        Region connectedRegion(*this);

        connectedRegion.startSlices = region.startSlices;
        connectedRegion.endSlices = region.endSlices;
        connectedRegion.signedNormals = region.signedNormals;

        std::queue<Coord> queue;

        queue.push(Coord(x, y, z));
        visited(x, y, z) = true;

        while (!queue.empty()) {
          const Coord current = queue.front();
          queue.pop();

          connectedRegion.voxels.push_back(
              voxelIds(current.x(), current.y(), current.z()));

          for (const auto &offset : neighbors) {
            const int nx2 = static_cast<int>(current.x()) + offset.x();
            const int ny2 = static_cast<int>(current.y()) + offset.y();
            const int nz2 = static_cast<int>(current.z()) + offset.z();

            if (nx2 < 0 || ny2 < 0 || nz2 < 0 || nx2 >= static_cast<int>(nx) ||
                ny2 >= static_cast<int>(ny) || nz2 >= static_cast<int>(nz)) {
              continue;
            }

            const size_t ux = static_cast<size_t>(nx2);
            const size_t uy = static_cast<size_t>(ny2);
            const size_t uz = static_cast<size_t>(nz2);

            if (!tensor(ux, uy, uz) || visited(ux, uy, uz))
              continue;

            visited(ux, uy, uz) = true;
            queue.push(Coord(ux, uy, uz));
          }
        }

        result.push_back(std::move(connectedRegion));
      }
    }
  }

  return result;
}

void VoxelMesh::splitByConnectedRegions() {
  std::vector<Region> res;
  for (RegionId rid = 0; rid < regions.size(); rid++) {
    auto ccs = getConnectedRegions(rid);
    // assert(ccs.size() == 1);
    res.insert(res.begin(), ccs.begin(), ccs.end());
  }
  regions = res;
}
