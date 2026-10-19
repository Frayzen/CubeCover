#ifndef FRAME_FIELD_HPP
#define FRAME_FIELD_HPP

#include <Eigen/Core>

Eigen::Matrix3i computePiST(Eigen::Matrix3d frameS, Eigen::Matrix3d frameT);

void fillRandomPermutations(std::vector<Eigen::Matrix3d> &frames);

#endif
