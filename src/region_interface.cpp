#include "region_interface.hpp"
#include <iostream>
#include <queue>

InterfacePosition
RegionInterface::getInterfacePosition(Axis cutDirection) const {
  assert(cutDirection != normal);
  Axis cutDirectionNextAxis = (cutDirection + 1) % DIMENSION;
  if (cutDirectionNextAxis == normal) {
    // LEFT RIGHT
    return sign == POSITIVE ? IP_LEFT : IP_RIGHT;
  } else {
    // TOP BOT
    return sign == POSITIVE ? IP_TOP : IP_BOT;
  }
}
static std::vector<Eigen::Vector2i>
getEdgeCells(InterfacePosition dir, Eigen::Index rows, Eigen::Index cols) {
  std::vector<Eigen::Vector2i> cells;

  switch (dir) {
  case IP_TOP:
    for (Eigen::Index col = 0; col < cols; ++col)
      cells.emplace_back(0, col);
    break;

  case IP_BOT:
    for (Eigen::Index col = 0; col < cols; ++col)
      cells.emplace_back(rows - 1, col);
    break;

  case IP_LEFT:
    for (Eigen::Index row = 0; row < rows; ++row)
      cells.emplace_back(row, 0);
    break;

  case IP_RIGHT:
    for (Eigen::Index row = 0; row < rows; ++row)
      cells.emplace_back(row, cols - 1);
    break;
  }

  return cells;
}

static bool
hasPathBetween(const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> &mat,
               const std::vector<Eigen::Vector2i> &starts,
               const std::vector<Eigen::Vector2i> &targets) {

  const Eigen::Index rows = mat.rows();
  const Eigen::Index cols = mat.cols();

  if (rows == 0 || cols == 0)
    return false;

  std::vector<std::vector<bool>> visited(rows, std::vector<bool>(cols, false));

  std::queue<Eigen::Vector2i> queue;

  for (const auto &p : starts) {
    if (mat(p.x(), p.y())) {
      queue.push(p);
      visited[p.x()][p.y()] = true;
    }
  }

  constexpr int dx[] = {-1, 1, 0, 0};
  constexpr int dy[] = {0, 0, -1, 1};

  while (!queue.empty()) {
    Eigen::Vector2i p = queue.front();
    queue.pop();

    for (const auto &target : targets) {
      if (p == target)
        return true;
    }

    for (int i = 0; i < 4; ++i) {
      Eigen::Index x = p.x() + dx[i];
      Eigen::Index y = p.y() + dy[i];

      if (x < 0 || x >= rows || y < 0 || y >= cols)
        continue;

      if (!visited[x][y] && mat(x, y)) {
        visited[x][y] = true;
        queue.emplace(x, y);
      }
    }
  }

  return false;
}

bool findPathInMatrix(
    const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> &mat,
    const InterfacePosition &from, const InterfacePosition &to) {

  assert(from != to);

  if (mat.rows() == 0 || mat.cols() == 0)
    return false;

  auto starts = getEdgeCells(from, mat.rows(), mat.cols());
  auto targets = getEdgeCells(to, mat.rows(), mat.cols());
  for (auto start : starts)
    std::cout << " START IS " << start << std::endl;
  for (auto end : targets)
    std::cout << " ENDS IS " << end << std::endl;

  return hasPathBetween(mat, starts, targets);
}

bool findPathInMatrixFrom(
    const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> &mat,
    const Eigen::Vector2i &from, const InterfacePosition &to) {

  if (mat.rows() == 0 || mat.cols() == 0)
    return false;

  auto starts = {from};
  auto targets = getEdgeCells(to, mat.rows(), mat.cols());

  return hasPathBetween(mat, starts, targets);
}
