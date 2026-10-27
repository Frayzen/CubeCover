#include "frame_field.hpp"
#include "constants.hpp"
#include "mesh_loader.hpp"

#include <Eigen/Eigen>
#include <array>
#include <cmath>

void fillRandomPermutations(std::vector<Eigen::Matrix3d> &frames,
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
    //   Vertex bary = mesh.tets[tet].getBarycenter(mesh.vertices);
    //   const double SQUARE_BOUND = 15;
    //   // const double SQUARE_BOUND = 0.15;
    //   if (bary.x() < SQUARE_BOUND && bary.y() < SQUARE_BOUND) {
    //     // case bot left
    //     frames[tet].col(0) << 1, 0, 0;
    //     frames[tet].col(1) << 0, 1, 0;
    //     frames[tet].col(2) << 0, 0, 1;
    //   } else if (bary.y() > SQUARE_BOUND && bary.x() < bary.y()) {
    //     // case top left
    //     double theta = M_PI / 3;
    //     frames[tet].col(0) << cos(theta), sin(theta), 0;
    //     frames[tet].col(1) << -sin(theta), cos(theta), 0;
    //     frames[tet].col(2) << 0, 0, 1;
    //   } else {
    //     // case bot right
    //     double theta = M_PI / 6;
    //     frames[tet].col(0) << cos(theta), sin(theta), 0;
    //     frames[tet].col(1) << -sin(theta), cos(theta), 0;
    //     frames[tet].col(2) << 0, 0, 1;
    //   }
    //   // auto frame = permutations[dist(rng)];
    frames[tet] = Eigen::Matrix3d::Identity();
  }
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
