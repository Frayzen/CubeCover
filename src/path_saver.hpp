#ifndef PATH_SAVER_HPP
#define PATH_SAVER_HPP

#include "region_divider.hpp"
#include <string>
#include <vector>

void savePaths(std::string filename, std::vector<Path> paths,
               std::vector<Vertex> vertices);

void savePaths(std::string filename, std::vector<std::vector<size_t>> paths,
               std::vector<Vertex> vertices);

#endif
