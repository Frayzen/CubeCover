#pragma once

#include <Eigen/Core>
#include <cstddef>
#include <sys/types.h>
using SliceId = size_t;
using NodeId = size_t;
using TetId = size_t;
using FaceId = size_t;
using EdgeId = size_t;
using RegionId = ssize_t;
const RegionId NO_REGION = -1;
using VoxelId = ssize_t;
const VoxelId NO_VOXEL = -1;
using RegionLayerId = ssize_t;
const RegionLayerId NO_LAYER = -1;
using Vertex = Eigen::Vector3d;
using Voxel = Eigen::Vector3i;
using ZoneId = size_t;
using LoopSegmentId = size_t;
