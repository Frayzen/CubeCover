#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cstddef>
#include <limits>

const double mind = std::numeric_limits<double>::min();
const double maxd = std::numeric_limits<double>::max();

const size_t DIMENSION = 3;
const size_t TET_NODE_COUNT = 4;
const size_t AMOUNT_SHAPE_FUNCTION = 4;

#endif
