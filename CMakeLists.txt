cmake_minimum_required (VERSION 2.8)

project (cloud_physics Fortran)

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
find_package (NETCDF REQUIRED COMPONENTS Fortran)
include_directories (${NETCDF_INCLUDE_DIRS})

include_directories ("${PROJECT_SOURCE_DIR}/src")

add_library (cloud_physics
    "${PROJECT_SOURCE_DIR}/src/constants.F90"
    "${PROJECT_SOURCE_DIR}/src/cloud_physics.F90"
    "${PROJECT_SOURCE_DIR}/src/advection_solver.F90"
    "${PROJECT_SOURCE_DIR}/src/lagrangian_solver.F90"
    "${PROJECT_SOURCE_DIR}/src/ui.F90"
)

add_executable (main "${PROJECT_SOURCE_DIR}/src/main.F90")
target_link_libraries (main cloud_physics ${NETCDF_LIBRARIES})
