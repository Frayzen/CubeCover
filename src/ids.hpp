#pragma once

#include <Eigen/Core>
#include <cstddef>
#include <sys/types.h>
using SliceId = int;

const SliceId MAX_SLICE = std::numeric_limits<SliceId>::max();
const SliceId MIN_SLICE = std::numeric_limits<SliceId>::min();

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
using PatchId = ssize_t;
const PatchId NO_PATCH = -1;

using GraphNodeId = ssize_t;
const GraphNodeId NO_NODE = -1;

using Vertex = Eigen::Vector3d;
using Voxel = Eigen::Vector3<SliceId>;
using ZoneId = size_t;
using LoopSegmentId = size_t;
