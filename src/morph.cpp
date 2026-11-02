#include "morph.hpp"

#include <vtkNew.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <stdexcept>

// alpha = 0 -> source mesh
// alpha = 1 -> target mesh
vtkSmartPointer<vtkUnstructuredGrid> morphMesh(const TetrahedralMesh &source,
                                               const TetrahedralMesh &target,
                                               double alpha) {
  // Check interpolation parameter
  alpha = std::clamp(alpha, 0.0, 1.0);

  // Check that both meshes have the same topology
  if (source.vertices.size() != target.vertices.size() ||
      source.tets.size() != target.tets.size()) {
    throw std::runtime_error("Meshes do not have the same topology.");
  }

  const size_t nbVertices = source.vertices.size();
  const size_t nbTets = source.tets.size();

  // Create VTK points
  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(nbVertices);

  // Interpolate vertex positions
  for (size_t i = 0; i < nbVertices; ++i) {
    const auto &p0 = source.vertices[i];
    const auto &p1 = target.vertices[i];

    double x = (1.0 - alpha) * p0.x() + alpha * p1.x();
    double y = (1.0 - alpha) * p0.y() + alpha * p1.y();
    double z = (1.0 - alpha) * p0.z() + alpha * p1.z();

    points->SetPoint(i, x, y, z);
  }

  // Create the VTK unstructured grid
  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(points);

  // Insert tetrahedral cells
  for (size_t i = 0; i < nbTets; ++i) {
    const auto &tet0 = source.tets[i];
    const auto &tet1 = target.tets[i];

    // Verify connectivity is identical
    for (int j = 0; j < 4; ++j) {
      if (tet0.nodes[j] != tet1.nodes[j]) {
        throw std::runtime_error(
            "Meshes have different tetrahedral connectivity.");
      }
    }

    vtkIdType ids[4];

    for (int j = 0; j < 4; ++j)
      ids[j] = static_cast<vtkIdType>(tet0.nodes[j]);

    grid->InsertNextCell(VTK_TETRA, 4, ids);
  }

  return grid;
}

void saveMorphAnimation(const TetrahedralMesh &source,
                        const TetrahedralMesh &target,
                        const std::string &outputName, int nbFrames) {
  for (int frame = 0; frame <= nbFrames; ++frame) {
    double alpha = static_cast<double>(frame) / nbFrames;

    auto morphed = morphMesh(source, target, alpha);

    std::ostringstream filename;
    filename << outputName << "_" << std::setw(4) << std::setfill('0') << frame
             << ".vtu";

    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(filename.str().c_str());
    writer->SetInputData(morphed);
    writer->SetDataModeToBinary();

    if (writer->Write() != 1)
      throw std::runtime_error("Failed to save frame");
  }

  // Create PVD file
  std::ofstream pvd(outputName + ".pvd");

  pvd << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"Collection\" version=\"0.1\" "
      << "byte_order=\"LittleEndian\">\n"
      << "<Collection>\n";

  for (int frame = 0; frame <= nbFrames; ++frame) {
    pvd << "<DataSet timestep=\"" << static_cast<double>(frame) / nbFrames
        << "\" group=\"\" part=\"0\" file=\"" << outputName << "_"
        << std::setw(4) << std::setfill('0') << frame << ".vtu\"/>\n";
  }

  pvd << "</Collection>\n"
      << "</VTKFile>\n";
}
