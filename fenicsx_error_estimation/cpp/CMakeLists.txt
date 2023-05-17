# Top level CMakeLists for fenicsx_error_estimation
cmake_minimum_required(VERSION 3.10)

# Set project name and version number
project(FENICSX_ERROR_ESTIMATION VERSION 0.5.0)

#------------------------------------------------------------------------------
# Set CMake options, see `cmake --help-policy CMP00xx`

cmake_policy(VERSION 3.10)
if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
if (POLICY CMP0075)
  cmake_policy(SET CMP0075 NEW)
endif()

# Use C++17
set(CMAKE_CXX_STANDARD 20)

# Require C++17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Do not enable compler-specific extensions
set(CMAKE_CXX_EXTENSIONS OFF)

#------------------------------------------------------------------------------
# configurable options for how we want to build

include(FeatureSummary)

option(BUILD_SHARED_LIBS "Build FENICSX_ERROR_ESTIMATION with shared libraries." ON)
add_feature_info(BUILD_SHARED_LIBS BUILD_SHARED_LIBS "Build FENICSX_ERROR_ESTIMATION with shared libraries.")

option(CMAKE_INSTALL_RPATH_USE_LINK_PATH "Add paths to linker search and installed rpath." ON)
add_feature_info(CMAKE_INSTALL_RPATH_USE_LINK_PATH CMAKE_INSTALL_RPATH_USE_LINK_PATH "Add paths to linker search and installed rpath.")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

# Find packages
find_package(DOLFINX REQUIRED)
find_package(Basix REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

feature_summary(WHAT ALL)

add_library(fenicsx_error_estimation "")

target_link_libraries(fenicsx_error_estimation PUBLIC Basix::basix)
target_link_libraries(fenicsx_error_estimation PUBLIC dolfinx)

# Installation of FENICSX_EESTM Library
include(GNUInstallDirs)

set(
  HEADERS_EESTM
  projected_local_solver.hpp
)

target_sources(
  fenicsx_error_estimation
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/projected_local_solver.cpp
)

install(FILES ${HEADERS_EESTM} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fenicsx_error_estimation COMPONENT Development)

# Set target include location (for build and installed)
target_include_directories(fenicsx_error_estimation PUBLIC
                           $<INSTALL_INTERFACE:include>
                           "$<BUILD_INTERFACE:${FENICSX_ERROR_ESTIMATION_SOURCE_DIR};${FENICSX_ERROR_ESTIMATION_SOURCE_DIR}/fenicsx_error_estimation>")

# Set target properties
set_target_properties(fenicsx_error_estimation PROPERTIES
  VERSION ${FENICSX_ERROR_ESTIMATION_VERSION}
  SOVERSION ${FENICSX_ERROR_ESTIMATION_VERSION_MAJOR}.${FENICSX_ERROR_ESTIMATION_VERSION_MINOR})

#------------------------------------------------------------------------------
# Install dolfinx library and header files
install(TARGETS fenicsx_error_estimation
  EXPORT FENICSX_ERROR_ESTIMATIONTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeExecutables
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
  )

# Generate DOLFINTargets.cmake
install(EXPORT FENICSX_ERROR_ESTIMATIONTargets DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/fenicsx_error_estimation)

include(CMakePackageConfigHelpers)

configure_package_config_file(${FENICSX_ERROR_ESTIMATION_SOURCE_DIR}/FENICSX_ERROR_ESTIMATIONConfig.cmake.in
${CMAKE_BINARY_DIR}/fenicsx_error_estimation/FENICSX_ERROR_ESTIMATIONConfig.cmake
INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/fenicsx_error_estimation)

# Install CMake helper files
install(
  FILES
  ${CMAKE_BINARY_DIR}/fenicsx_error_estimation/FENICSX_ERROR_ESTIMATIONConfig.cmake
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/fenicsx_error_estimation
  COMPONENT Development)