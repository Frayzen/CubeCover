#ifndef MESH_SAVER_HPP
#define MESH_SAVER_HPP

#include "mesh_loader.hpp"
#include "vertex.hpp"

#include "region_divider.hpp"
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

enum class ElementType { Node, Edge, Cell };

std::vector<Edge> buildEdges(const std::vector<Tetrahedron> &tets);

class MeshSaver {
private:
  vtkSmartPointer<vtkUnstructuredGrid> grid;
  vtkSmartPointer<vtkPolyData> edgeGrid;

  std::string path;

  size_t nbNodes;
  size_t nbEdges;
  size_t nbCells;

  template <typename T> static constexpr int nbComponents() {
    if constexpr (std::is_same_v<T, bool>) {
      return 1;
    } else if constexpr (std::is_arithmetic_v<T>) {
      return 1;
    } else if constexpr (requires {
                           T::RowsAtCompileTime;
                           T::ColsAtCompileTime;
                         }) {
      return T::RowsAtCompileTime * T::ColsAtCompileTime;
    } else {
      static_assert(sizeof(T) == 0, "Unsupported field type");
    }
  }

  template <typename T> static double component(const T &value, size_t i) {
    if constexpr (std::is_arithmetic_v<T>) {
      return static_cast<double>(value);
    } else {
      constexpr size_t rows = T::RowsAtCompileTime;
      return static_cast<double>(value(i % rows, i / rows));
    }
  }

  template <typename T>
  void addFieldToArray(vtkDataSetAttributes *attributes,
                       const std::string &name, const std::vector<T> &values) {
    auto array = vtkSmartPointer<vtkDoubleArray>::New();

    constexpr int nComp = nbComponents<T>();

    array->SetName(name.c_str());
    array->SetNumberOfComponents(nComp);
    array->SetNumberOfTuples(values.size());

    for (size_t i = 0; i < values.size(); ++i) {
      double tuple[nComp];

      for (int j = 0; j < nComp; ++j) {
        tuple[j] = component(values[i], j);
      }

      array->SetTuple(i, tuple);
    }

    attributes->AddArray(array);
  }

public:
  MeshSaver(const std::vector<Vertex> &vertices, const std::vector<Edge> &edges,
            const std::vector<Tetrahedron> &tets, const std::string &path);

  MeshSaver(const std::vector<Vertex> &vertices,
            const std::vector<Tetrahedron> &tets, const std::string &path)
      : MeshSaver(vertices, buildEdges(tets), tets, path) {}

  template <ElementType E, typename T>
  void addField(const std::string &name, const std::vector<T> &values) {
    vtkDataSetAttributes *attributes = nullptr;
    size_t expectedSize = 0;

    if constexpr (E == ElementType::Node) {
      attributes = grid->GetPointData();
      expectedSize = nbNodes;
    } else if constexpr (E == ElementType::Cell) {
      attributes = grid->GetCellData();
      expectedSize = nbCells;
    } else if constexpr (E == ElementType::Edge) {
      attributes = edgeGrid->GetCellData();
      expectedSize = nbEdges;
    }

    if (values.size() != expectedSize) {
      throw std::invalid_argument(
          "Field size does not match the number of elements.");
    }

    addFieldToArray(attributes, name, values);
  }

  void saveMesh() const;

  void saveEdges() const;
};

void savePaths(std::string filename, std::vector<Path> paths,
               std::vector<Vertex> vertices);

void savePaths(std::string filename, std::vector<std::vector<size_t>> paths,
               std::vector<Vertex> vertices);

#endif
