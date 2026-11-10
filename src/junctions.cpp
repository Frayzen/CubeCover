#include "junctions.hpp"
#include "old/signed_axis.hpp"
#include "voxel_mesh.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_set>
#include <vector>

struct FacePlane {
  Vertex normal;
  double offset;
};

const double relativeTolerance = 0.04;
const double normalTolerance = 1e-1;
const double planeTolerance = 1e-0;
const size_t minComponentSize = 3;

// -----------------------------------------------------------------------------
// Surface patch construction
// -----------------------------------------------------------------------------

static FacePlane computeFacePlane(FaceId faceId, TetrahedralMesh &mesh) {
  const Face &face = mesh.faces[faceId];

  const Vertex &a = mesh.vertices[face[0]];
  const Vertex &b = mesh.vertices[face[1]];
  const Vertex &c = mesh.vertices[face[2]];

  Vertex normal = (b - a).cross(c - a);
  const double length = normal.norm();

  assert(length > 1e-14);

  normal /= length;

  return {normal, normal.dot(a)};
}

static bool areCoplanar(const FacePlane &a, const FacePlane &b,
                        double normalTolerance, double planeTolerance) {
  const double alignment = a.normal.dot(b.normal);

  if (std::abs(alignment) < 1.0 - normalTolerance)
    return false;

  const double sign = alignment >= 0.0 ? 1.0 : -1.0;

  return std::abs(a.offset - sign * b.offset) <= planeTolerance;
}

SurfacePatches buildSurfacePatches(TetrahedralMesh &mesh) {
  SurfacePatches result = {mesh};

  result.faceToPatch.assign(mesh.nbBoundaryFace, -1);

  std::vector<FacePlane> planes(mesh.nbBoundaryFace);

  for (FaceId faceId = 0; faceId < mesh.nbBoundaryFace; ++faceId) {
    planes[faceId] = computeFacePlane(faceId, mesh);
  }

  for (FaceId seed = 0; seed < mesh.nbBoundaryFace; ++seed) {
    if (result.faceToPatch[seed] != -1)
      continue;

    const int patchId = static_cast<int>(result.patches.size());

    result.patches.emplace_back();

    std::queue<FaceId> queue;
    queue.push(seed);
    result.faceToPatch[seed] = patchId;

    while (!queue.empty()) {
      const FaceId current = queue.front();
      queue.pop();

      result.patches.back().push_back(current);

      for (size_t edgeId : mesh.faceToEdges[current]) {
        const auto &adjacentFaces = mesh.edgeToBoundaryFaces[edgeId];

        for (ssize_t adjacentFace : adjacentFaces) {
          if (adjacentFace < 0 ||
              adjacentFace >= static_cast<ssize_t>(mesh.nbBoundaryFace))
            continue;

          const FaceId neighbor = static_cast<FaceId>(adjacentFace);

          if (result.faceToPatch[neighbor] != -1)
            continue;

          if (!areCoplanar(planes[current], planes[neighbor], normalTolerance,
                           planeTolerance))
            continue;

          result.faceToPatch[neighbor] = patchId;
          queue.push(neighbor);
        }
      }
    }
  }

  return result;
}

// -----------------------------------------------------------------------------
// Point-to-triangle closest point
// -----------------------------------------------------------------------------

static Vertex closestPointOnTriangle(const Vertex &p, const Vertex &a,
                                     const Vertex &b, const Vertex &c) {
  const Vertex ab = b - a;
  const Vertex ac = c - a;
  const Vertex ap = p - a;

  const double d1 = ab.dot(ap);
  const double d2 = ac.dot(ap);

  if (d1 <= 0.0 && d2 <= 0.0)
    return a;

  const Vertex bp = p - b;
  const double d3 = ab.dot(bp);
  const double d4 = ac.dot(bp);

  if (d3 >= 0.0 && d4 <= d3)
    return b;

  const double vc = d1 * d4 - d3 * d2;

  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    const double v = d1 / (d1 - d3);
    return a + v * ab;
  }

  const Vertex cp = p - c;
  const double d5 = ab.dot(cp);
  const double d6 = ac.dot(cp);

  if (d6 >= 0.0 && d5 <= d6)
    return c;

  const double vb = d5 * d2 - d1 * d6;

  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    const double w = d2 / (d2 - d6);
    return a + w * ac;
  }

  const double va = d3 * d6 - d5 * d4;

  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));

    return b + w * (c - b);
  }

  const double denominator = 1.0 / (va + vb + vc);
  const double v = vb * denominator;
  const double w = vc * denominator;

  return a + v * ab + w * ac;
}

static double pointToPatchDistanceSquared(const Vertex &point, int patchId,
                                          const SurfacePatches &surface) {
  auto &mesh = surface.mesh;
  double minDistanceSquared = std::numeric_limits<double>::infinity();

  for (FaceId faceId : surface.patches[patchId]) {
    const Face &face = mesh.faces[faceId];

    const Vertex closest =
        closestPointOnTriangle(point, mesh.vertices[face[0]],
                               mesh.vertices[face[1]], mesh.vertices[face[2]]);

    const double distanceSquared = (point - closest).squaredNorm();

    minDistanceSquared = std::min(minDistanceSquared, distanceSquared);
  }

  return minDistanceSquared;
}

// -----------------------------------------------------------------------------
// Main function
// -----------------------------------------------------------------------------

std::vector<Voxel> computeMedialAxisJunctions(const VoxelMesh &mesh,
                                              const SurfacePatches &surface) {

  Eigen::DSizes<Eigen::Index, 3> dim = mesh.voxels.dimensions();

  for (Axis d = 0; d < DIMENSION; d++)
    dim[d]++;

  Eigen::Tensor<PatchId, DIMENSION> nodeClosestPatch(dim);

  // Compute closest surface patch for every grid node
  for (size_t i = 0; i <= mesh.voxels.dimension(0); i++) {
    for (size_t j = 0; j <= mesh.voxels.dimension(1); j++) {
      for (size_t k = 0; k <= mesh.voxels.dimension(2); k++) {

        Vertex v = mesh.origin + mesh.voxelSize * Vertex(i, j, k);

        double minDist = INFINITY;
        PatchId id = NO_PATCH;

        for (size_t p = 0; p < surface.patches.size(); p++) {
          double dist = pointToPatchDistanceSquared(v, p, surface);

          if (dist < minDist) {
            minDist = dist;
            id = static_cast<PatchId>(p);
          }
        }

        nodeClosestPatch(i, j, k) = id;
      }
    }
  }

  const std::vector<Voxel> offsets = {
      {-1, 1, 1},  {1, 1, 1},  {-1, -1, 1},  {1, -1, 1},
      {-1, 1, -1}, {1, 1, -1}, {-1, -1, -1}, {1, -1, -1},
  };

  std::vector<Voxel> junctions;
  for (int i = 0; i < mesh.voxels.dimension(0); i++) {
    for (int j = 0; j < mesh.voxels.dimension(1); j++) {
      for (int k = 0; k < mesh.voxels.dimension(2); k++) {
        if (mesh.voxels(i, j, k) == NO_VOXEL)
          continue;
        std::unordered_set<PatchId> closestPatches;
        for (const Voxel &offset : offsets) {
          PatchId id = nodeClosestPatch(
              i + (offset.x() > 0), j + (offset.y() > 0), k + (offset.z() > 0));
          if (id != NO_PATCH)
            closestPatches.insert(id);
        }
        if (closestPatches.size() >= 3)
          junctions.emplace_back(i, j, k);
      }
    }
  }

  return junctions;
}

std::vector<Vertex>
computeJunctionClusters(TetrahedralMesh &mesh,
                        const std::vector<TetId> &junctionTetIds) {

  std::vector<Vertex> clusterBarycenters;

  if (junctionTetIds.empty())
    return clusterBarycenters;

  // Mark which tetrahedra belong to the input list.
  std::vector<bool> isJunctionTet(mesh.tets.size(), false);

  for (TetId tetId : junctionTetIds) {
    assert(tetId < mesh.tets.size());
    isJunctionTet[tetId] = true;
  }

  // For each mesh vertex, store the candidate tetrahedra containing it.
  std::vector<std::vector<TetId>> nodeToJunctionTets(mesh.vertices.size());

  for (TetId tetId : junctionTetIds) {
    const Tetrahedron &tet = mesh.tets[tetId];

    for (size_t localId = 0; localId < NB_NODE_PER_TET; ++localId) {
      nodeToJunctionTets[tet[localId]].push_back(tetId);
    }
  }

  // Find connected components.
  std::vector<bool> visited(mesh.tets.size(), false);

  for (TetId seed : junctionTetIds) {
    if (visited[seed])
      continue;

    std::queue<TetId> queue;
    queue.push(seed);
    visited[seed] = true;

    Vertex barycenter = Vertex::Zero();
    size_t nbTets = 0;

    while (!queue.empty()) {
      const TetId current = queue.front();
      queue.pop();

      // Accumulate the barycenter of the current tetrahedron.
      barycenter += mesh.tets[current].getBarycenter(mesh.vertices);
      ++nbTets;

      // Find all candidate tetrahedra sharing one of its vertices.
      const Tetrahedron &tet = mesh.tets[current];

      for (size_t localId = 0; localId < NB_NODE_PER_TET; ++localId) {
        const NodeId nodeId = tet[localId];

        for (TetId neighbor : nodeToJunctionTets[nodeId]) {
          if (visited[neighbor])
            continue;

          visited[neighbor] = true;
          queue.push(neighbor);
        }
      }
    }

    barycenter /= static_cast<double>(nbTets);
    clusterBarycenters.push_back(barycenter);
  }

  return clusterBarycenters;
}
