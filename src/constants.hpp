#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cstddef>
#include <limits>
#include <sys/types.h>

const double mind = std::numeric_limits<double>::min();
const double maxd = std::numeric_limits<double>::max();

const size_t DIMENSION = 3;
const size_t NB_NODE_PER_FACE = 3;
const size_t NB_NODE_PER_TET = 4;
const size_t NB_NODE_PER_HEX = 8;
const size_t NB_UNKNOWN_PER_TET = NB_NODE_PER_TET * DIMENSION;

const ssize_t UNDEFINED_ELEMENT = -1;

const int INT_NAN = std::numeric_limits<int>::max();

const double MIN_QUALITY_INCREASE_PER_LOOP = 0.1;

#endif
