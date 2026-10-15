#include "mesh_loader.hpp"
#include "constants.hpp"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

static vtkSmartPointer<vtkUnstructuredGrid> readGrid(const std::string &path) {
  const std::filesystem::path filePath(path);
  const std::string extension = filePath.extension().string();

  if (extension == ".vtu" || extension == ".VTU") {
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(path.c_str());
    reader->Update();

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->DeepCopy(reader->GetOutput());
    return grid;
  }

  if (extension == ".vtk" || extension == ".VTK") {
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(path.c_str());
    reader->Update();

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->DeepCopy(reader->GetOutput());
    return grid;
  }

  std::cerr << "Error: unsupported file extension '" << extension
            << "'. Supported inputs are .vtu and legacy .vtk "
               "unstructured-grid files.\n";
  return nullptr;
}

#include <algorithm>
#include <map>

bool loadTetrahedralMesh(const std::string &path, TetrahedralMesh &output) {
  auto grid = readGrid(path);
  if (!grid) {
    return false;
  }

  if (grid->GetNumberOfPoints() == 0) {
    std::cerr << "Error: mesh contains no points or could not be read.\n";
    return false;
  }

  if (grid->GetNumberOfCells() == 0) {
    std::cerr << "Error: mesh contains no cells or could not be read.\n";
    return false;
  }

  // Validate all cells first
  for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
    const size_t cellType = grid->GetCellType(cellId);
    if (cellType != VTK_TETRA) {
      std::cerr << "Error: cell " << cellId << " has VTK cell type " << cellType
                << "; expected VTK_TETRA (" << VTK_TETRA << ").\n";
      return false;
    }

    vtkCell *cell = grid->GetCell(cellId);
    if (cell == nullptr || cell->GetNumberOfPoints() != 4) {
      std::cerr << "Error: cell " << cellId
                << " is not a valid four-node tetrahedron.\n";
      return false;
    }
  }

  TetrahedralMesh candidate;
  candidate.amountInternalFace = 0;
  vtkPoints *points = grid->GetPoints();
  candidate.vertices.reserve(
      static_cast<std::size_t>(grid->GetNumberOfPoints()));

  for (vtkIdType pointId = 0; pointId < grid->GetNumberOfPoints(); ++pointId) {
    double p[3];
    points->GetPoint(pointId, p);
    candidate.vertices.push_back(Vertex{p[0], p[1], p[2]});
  }

  candidate.tets.reserve(static_cast<std::size_t>(grid->GetNumberOfCells()));

  // Map to keep track of unique faces: sorted face vertex IDs -> list of
  // adjacent cell IDs
  std::map<std::array<vtkIdType, 3>, std::array<ssize_t, 2>> faceMap;

  for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
    vtkCell *cell = grid->GetCell(cellId);
    Tetrahedron tet{};

    for (size_t localId = 0; localId < 4; ++localId) {
      const vtkIdType pointId = cell->GetPointId(localId);
      if (pointId < 0 || pointId >= grid->GetNumberOfPoints()) {
        std::cerr << "Error: tetrahedron " << cellId
                  << " references an invalid posize_t index.\n";
        return false;
      }
      tet[static_cast<std::size_t>(localId)] = pointId;
    }

    candidate.tets.push_back(tet);

    // Extract the 4 faces of the tetrahedron and register cell adjacency
    for (size_t faceId = 0; faceId < cell->GetNumberOfFaces(); ++faceId) {
      vtkCell *face = cell->GetFace(faceId);
      if (face != nullptr && face->GetNumberOfPoints() == 3) {
        std::array<vtkIdType, 3> f_keys = {
            face->GetPointId(0), face->GetPointId(1), face->GetPointId(2)};

        // the key is independent of traversal order
        std::sort(f_keys.begin(), f_keys.end());

        if (faceMap.find(f_keys) == faceMap.end()) {
          faceMap[f_keys][0] = cellId;
          faceMap[f_keys][1] = UNDEFINED_ADJACENT_CELL;
        } else {
          candidate.amountInternalFace++;
          faceMap[f_keys][1] = cellId;
        }
      }
    }
  }

  // Populate candidate faces from the aggregated map
  candidate.faces.reserve(faceMap.size());
  for (const auto &[f_keys, adjacentCells] : faceMap) {
    candidate.faces.push_back(Face{f_keys, adjacentCells});
  }

  output = std::move(candidate);
  return true;
}
