cmake_minimum_required(VERSION 3.16)
project(VtkTetReader LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    IOLegacy
    IOXML
)

find_package(Eigen3 REQUIRED)

file(GLOB_RECURSE SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")
message(DEBUG "SOURCES: ${SOURCE_FILES}")

add_executable(CubeCover ${SOURCE_FILES})
target_include_directories(CubeCover PRIVATE ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(CubeCover PRIVATE Eigen3::Eigen)
target_link_libraries(CubeCover PRIVATE
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::IOLegacy
    VTK::IOXML
)

# find_package(Kokkos 4.2 REQUIRED CONFIG)

vtk_module_autoinit(
    TARGETS CubeCover
    MODULES
        VTK::CommonCore
        VTK::CommonDataModel
        VTK::IOLegacy
        VTK::IOXML
)

