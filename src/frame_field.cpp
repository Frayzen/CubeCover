#include "frame_field.hpp"
#include "constants.hpp"
#include "mesh_loader.hpp"
#include "trilaterate.hpp"

#include <Eigen/Eigen>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>

void fillFrameField(std::vector<Eigen::Matrix3d> &frames,
                    TetrahedralMesh &mesh) {

  std::array<Eigen::Matrix3d, 24> permutations;
  int count = 0;

  for (int p0 = 0; p0 < 3; ++p0) {
    for (int p1 = 0; p1 < 3; ++p1) {
      if (p1 == p0)
        continue;

      for (int p2 = 0; p2 < 3; ++p2) {
        if (p2 == p0 || p2 == p1)
          continue;

        for (int s0 : {-1, 1}) {
          for (int s1 : {-1, 1}) {
            for (int s2 : {-1, 1}) {
              Eigen::Matrix3d P = Eigen::Matrix3d::Zero();
              P(p0, 0) = s0;
              P(p1, 1) = s1;
              P(p2, 2) = s2;

              if (std::lround(P.determinant()) == 1)
                permutations[count++] = P;
            }
          }
        }
      }
    }
  }

  // std::mt19937 rng(std::random_device{}());
  // std::mt19937 rng(42);
  // std::uniform_int_distribution<int> dist(0, 23);

  for (size_t tet = 0; tet < mesh.tets.size(); tet++) {
    // auto facesId = mesh.tetToFaces[tet];
    // bool boundaryTet = false;
    // Vertex normalAcc = Vertex::Zero();
    // for (size_t faceId : facesId) {
    //   Face &face = mesh.faces[faceId];
    //   if (face.adjacentCells[1] == UNDEFINED_ELEMENT) {
    //     bool boundaryTet = true;
    //     auto vIds = face.vertices;
    //     Vertex nodes[3];
    //     for (size_t i = 0; i < NB_NODE_PER_FACE; i++)
    //       nodes[i] = mesh.vertices[vIds[i]];
    //     auto tetIds = mesh.tets[tet];
    //     size_t tipLocalId = NB_NODE_PER_TET;
    //     for (tipLocalId = 0; tipLocalId < NB_NODE_PER_TET; tipLocalId++) {
    //       size_t tipGlobalId = tetIds[tipLocalId];
    //       if (std::find(vIds.begin(), vIds.end(), tipGlobalId) == vIds.end())
    //         break;
    //     }
    //     assert(tipLocalId < NB_NODE_PER_TET);
    //     Vertex tip = mesh.vertices[tetIds[tipLocalId]];
    //     Vertex n = getNormal(nodes);
    //     if ((nodes[0] - tip).dot(n) < 0)
    //       n *= -1;
    //     normalAcc += n;
    //   }
    // }
    // if (boundaryTet) {
    //   normalAcc.normalize();

    //   size_t coordAlignBest;
    //   double signBest;
    //   double measureBest = 0;

    //   size_t coordAlignWorst;
    //   double signWorst;
    //   double measureWorst = std::numeric_limits<double>::max();

    //   for (size_t coord = 0; coord < DIMENSION; coord++) {
    //     Vertex unit = Vertex::Zero();
    //     unit[coord] = 1;
    //     double dotProduct = normalAcc.dot(unit);
    //     double sign = dotProduct > 0 ? 1 : -1;
    //     double cur = std::abs(dotProduct);
    //     if (measureBest + 10e-3 < cur) {
    //       measureBest = cur;
    //       coordAlignBest = coord;
    //       signBest = sign;
    //     }
    //     if (measureWorst - 10e-3 > cur) {
    //       measureWorst = cur;
    //       coordAlignWorst = coord;
    //       signWorst = sign;
    //     }
    //   }
    //   assert(measureWorst < 10e-3); // temporary, just for prism

    //   Eigen::Matrix3d frame;
    //   frame.col(coordAlignBest) = normalAcc;
    //   frame.col(coordAlignWorst) = Vertex::Zero();
    //   frame.col(coordAlignWorst)[coordAlignWorst] = 1;
    //   size_t lastCoord = 2;
    //   if ((coordAlignWorst + coordAlignBest) == 3)
    //     lastCoord = 0;
    //   if ((coordAlignWorst + coordAlignBest) == 2)
    //     lastCoord = 1;
    //   frame.col(lastCoord) =
    //   normalAcc.cross(frames[tet].col(coordAlignWorst)); frames[tet] +=
    //   frame; for (size_t coord = 0; coord < DIMENSION; coord++)
    //     frames[tet].col(coord).normalize();
    // } else {
    Vertex bary = Vertex::Zero();
    for (size_t i = 0; i < NB_NODE_PER_TET; i++) {
      bary += mesh.vertices[mesh.tets[tet][i]];
    }
    bary /= NB_NODE_PER_TET;

    if (bary.x() < 0.15 && bary.y() < 0.15) {
      // case bot left
      frames[tet].col(0) << 1, 0, 0;
      frames[tet].col(1) << 0, 1, 0;
      frames[tet].col(2) << 0, 0, 1;
    } else if (bary.y() > 0.15 && bary.x() < bary.y()) {
      // case top left
      double theta = M_PI / 3;
      frames[tet].col(0) << cos(theta), sin(theta), 0;
      frames[tet].col(1) << -sin(theta), cos(theta), 0;
      frames[tet].col(2) << 0, 0, 1;
    } else {
      // case bot right
      double theta = M_PI / 6;
      frames[tet].col(0) << cos(theta), sin(theta), 0;
      frames[tet].col(1) << -sin(theta), cos(theta), 0;
      frames[tet].col(2) << 0, 0, 1;
    }
  }
  // }
}

Eigen::Matrix3i computePiST(Eigen::Matrix3d frameS, Eigen::Matrix3d frameT) {
  // std::cout << "FRAME S " << std::endl << frameS << std::endl;
  // std::cout << "FRAME T " << std::endl << frameS << std::endl;
  double bestError = std::numeric_limits<double>::infinity();
  Eigen::Matrix3i bestPi = Eigen::Matrix3i::Identity();

  for (int p0 = 0; p0 < DIMENSION; ++p0) {
    for (int p1 = 0; p1 < DIMENSION; ++p1) {
      if (p1 == p0)
        continue;
      for (int p2 = 0; p2 < DIMENSION; ++p2) {
        if (p2 == p0 || p2 == p1)
          continue;

        Eigen::Matrix3i P = Eigen::Matrix3i::Zero();
        P(p0, 0) = 1;
        P(p1, 1) = 1;
        P(p2, 2) = 1;

        // Enumerate all sign combinations.
        for (int s0 : {-1, 1}) {
          for (int s1 : {-1, 1}) {
            for (int s2 : {-1, 1}) {
              Eigen::Matrix3i Pi = P;
              Pi.col(0) *= s0;
              Pi.col(1) *= s1;
              Pi.col(2) *= s2;

              // Keep only orientation-preserving matchings.
              if (Pi.determinant() != 1)
                continue;

              const Eigen::Matrix3d diff = frameS - frameT * Pi.cast<double>();

              const double error = diff.squaredNorm();

              if (error < bestError) {
                bestError = error;
                bestPi = Pi;
              }
            }
          }
        }
      }
    }
  }
  // std::cout << "BEST PI" << std::endl << bestPi << std::endl;
  return bestPi;
}
