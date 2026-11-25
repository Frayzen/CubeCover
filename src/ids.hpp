#pragma once

#include <Eigen/Core>
#include <cstddef>
#include <sys/types.h>
using SliceId = int;

const SliceId MAX_SLICE = std::numeric_limits<SliceId>::max();
const SliceId MIN_SLICE = std::numeric_limits<SliceId>::min();
const SliceId NO_SLICE = MIN_SLICE;

using SurfaceId = ssize_t;
const SurfaceId NO_SURFACE = -1;

using PathId = size_t;
using QuadId = size_t;
using NodeId = size_t;
using TetId = size_t;
using FaceId = size_t;
using EdgeId = ssize_t;
const EdgeId NO_EDGE = -1;

using PrimalRegionId = ssize_t;
const PrimalRegionId NO_PRIMAL_REGION = -1;

using DualRegionId = ssize_t;
const DualRegionId NO_DUAL_REGION = -1;
using VoxelId = ssize_t;
const VoxelId NO_VOXEL = -1;
using PatchId = ssize_t;
const PatchId NO_PATCH = -1;

using GraphNodeId = ssize_t;
const GraphNodeId NO_NODE = -1;

using Vertex = Eigen::Vector3d;
using Voxel = Eigen::Vector3<SliceId>;
using ZoneId = size_t;
using LoopSegmentId = size_t;

using FlagMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
