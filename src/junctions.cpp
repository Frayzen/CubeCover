#include "junctions.hpp"
#include "old/signed_axis.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <queue>
#include <vector>

struct FacePlane {
  Vertex normal;
  double offset;
};

const double relativeTolerance = 0.04;
const double normalTolerance = 1e-1;
const double planeTolerance = 1e-0;
const size_t minComponentSize = 3;

using NodePatchDistances = std::vector<std::vector<double>>;

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
  SurfacePatches result;

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
                                          const TetrahedralMesh &mesh,
                                          const SurfacePatches &surface) {
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
// Compute distances from every mesh vertex to every surface patch
// -----------------------------------------------------------------------------

static NodePatchDistances
computeNodePatchDistances(const TetrahedralMesh &mesh,
                          const SurfacePatches &surface) {
  NodePatchDistances distances(
      mesh.vertices.size(),
      std::vector<double>(surface.patches.size(),
                          std::numeric_limits<double>::infinity()));

  for (NodeId nodeId = 0; nodeId < mesh.vertices.size(); ++nodeId) {
    const Vertex &point = mesh.vertices[nodeId];

    for (int patchId = 0; patchId < static_cast<int>(surface.patches.size());
         ++patchId) {
      distances[nodeId][patchId] =
          std::sqrt(pointToPatchDistanceSquared(point, patchId, mesh, surface));
    }
  }

  return distances;
}

// -----------------------------------------------------------------------------
// Remove small connected components of candidate tetrahedra
// -----------------------------------------------------------------------------

static std::vector<TetId>
removeSmallComponents(const TetrahedralMesh &mesh,
                      const std::vector<bool> &candidates,
                      size_t minComponentSize) {
  std::vector<bool> visited(mesh.tets.size(), false);
  std::vector<TetId> result;

  for (TetId seed = 0; seed < mesh.tets.size(); ++seed) {
    if (!candidates[seed] || visited[seed])
      continue;

    std::vector<TetId> component;
    std::queue<TetId> queue;

    queue.push(seed);
    visited[seed] = true;

    while (!queue.empty()) {
      const TetId current = queue.front();
      queue.pop();

      component.push_back(current);

      // Two tetrahedra are adjacent if they share a face.
      for (ssize_t faceId : mesh.tetToFaces[current]) {
        if (faceId < 0)
          continue;

        const auto &adjacent = mesh.faces[faceId].adjacentCells;

        for (ssize_t neighbor : adjacent) {
          if (neighbor < 0 ||
              neighbor >= static_cast<ssize_t>(mesh.tets.size()))
            continue;

          const TetId neighborId = static_cast<TetId>(neighbor);

          if (!candidates[neighborId] || visited[neighborId])
            continue;

          visited[neighborId] = true;
          queue.push(neighborId);
        }
      }
    }

    if (component.size() >= minComponentSize) {
      result.insert(result.end(), component.begin(), component.end());
    }
  }

  return result;
}

// -----------------------------------------------------------------------------
// Main function
// -----------------------------------------------------------------------------

std::vector<TetId> computeMedialAxisJunctions(const TetrahedralMesh &mesh,
                                              const SurfacePatches &surface) {
  assert(mesh.nbBoundaryFace > 0);
  assert(relativeTolerance >= 0.0);

  // 2. Compute vertex-to-patch distances.
  const NodePatchDistances distances = computeNodePatchDistances(mesh, surface);

  // 3. Detect candidate tetrahedra.
  std::vector<bool> candidates(mesh.tets.size(), false);

  for (TetId tetId = mesh.nbBoundaryTet; tetId < mesh.tets.size(); ++tetId) {
    const Tetrahedron &tet = mesh.tets[tetId];

    std::vector<int> patchSupport(surface.patches.size(), 0);

    for (size_t localId = 0; localId < NB_NODE_PER_TET; ++localId) {
      const NodeId nodeId = tet[localId];
      const auto &nodeDistances = distances[nodeId];

      const double minDistance =
          *std::min_element(nodeDistances.begin(), nodeDistances.end());

      if (!std::isfinite(minDistance))
        continue;

      for (size_t patchId = 0; patchId < nodeDistances.size(); ++patchId) {
        if (nodeDistances[patchId] <= (1.0 + relativeTolerance) * minDistance) {
          ++patchSupport[patchId];
        }
      }
    }

    // A patch must be competitive at least two tet vertices.
    size_t supportedPatches = 0;

    for (int support : patchSupport) {
      if (support >= 2)
        ++supportedPatches;
    }

    candidates[tetId] = supportedPatches >= 3;
  }

  // 4. Remove isolated or very small candidate components.
  return removeSmallComponents(mesh, candidates, minComponentSize);
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
