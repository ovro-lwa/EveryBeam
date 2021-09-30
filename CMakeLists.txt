# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

#------------------------------------------------------------------------------
# Top level CMakeLists.txt file for EveryBeam
cmake_minimum_required(VERSION 3.7)

#------------------------------------------------------------------------------
# Set version name and project number
set(EVERYBEAM_VERSION 0.2.2)
if(EVERYBEAM_VERSION MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+)")
  set(EVERYBEAM_VERSION_MAJOR "${CMAKE_MATCH_1}")
  set(EVERYBEAM_VERSION_MINOR "${CMAKE_MATCH_2}")
  set(EVERYBEAM_VERSION_PATCH "${CMAKE_MATCH_3}")
else()
  message(
    FATAL_ERROR "Failed to parse EVERYBEAM_VERSION='${EVERYBEAM_VERSION}'")
endif()

project(EveryBeam VERSION ${EVERYBEAM_VERSION})

option(BUILD_WITH_PYTHON "Build python bindings" OFF)
option(BUILD_APT_PACKAGES "Build apt package" OFF)
option(BUILD_TESTING "Build tests" OFF)

string(TOLOWER ${CMAKE_PROJECT_NAME} projectname)

# Set the path to CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

# Configure directory for data files
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -DNDEBUG")
string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
set(CMAKE_INSTALL_DATA_DIR "${CMAKE_INSTALL_PREFIX}/share/${projectname}")
message("Storing data in: " ${CMAKE_INSTALL_DATA_DIR})

# Find and include git submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  # Update submodules as needed
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --checkout
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(
        FATAL_ERROR
          "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules"
      )
    endif()
  endif()
endif()

# Include aocommon/eigen3/pybind11 headers
include_directories("${CMAKE_SOURCE_DIR}/external/aocommon/include/")
include_directories("${CMAKE_SOURCE_DIR}/external/eigen/")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/pybind11")
include_directories(${pybind11_INCLUDE_DIR})

# Include schaapcommon
add_subdirectory("${CMAKE_SOURCE_DIR}/external/schaapcommon")
include_directories("${CMAKE_SOURCE_DIR}/external/schaapcommon/include")

# Find and include HDF5
find_package(
  HDF5
  COMPONENTS C CXX
  REQUIRED)
add_definitions(${HDF5_DEFINITIONS} -DH5_USE_110_API)
include_directories(${HDF5_INCLUDE_DIR})

# Find and include Casacore
set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)
include_directories(${CASACORE_INCLUDE_DIR})

# Find and include OpenMP
find_package(OpenMP REQUIRED)

# Find and include Boost headers (boost::math required for MWA beam)
# Filesystem component required to search for LOBES coefficient files
find_package(
  Boost
  COMPONENTS filesystem system
  REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Find and include FFTW3 float libraries
find_package(
  FFTW3
  COMPONENTS single
  REQUIRED)

# Find CFITSIO, needed for everybeam::aterms only
find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIR})

#------------------------------------------------------------------------------
# Set CMake and compiler options

# wget needs to be installed in order to download coefficient files on the fly
include(FindWget)
if(NOT WGET_FOUND)
  message(FATAL_ERROR "wget not found. Install wget on your local system")
endif()

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

# Set compile options
add_compile_options(-std=c++11 "${OpenMP_CXX_FLAGS}" -Wall)
string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")

if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
  add_compile_options(-DNDEBUG)
endif()

# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
     "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

#------------------------------------------------------------------------------
# Add source
add_subdirectory(cpp)

#------------------------------------------------------------------------------
# Add tests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  include(CTest)

  add_subdirectory(cpp/test)

  # TODO: compiling the demos should probably be a different cmake project
  # in which we use find_package(EveryBeam)
  add_subdirectory(demo)
endif()

#------------------------------------------------------------------------------
# Generate config.h and version.h headers
configure_file(${CMAKE_SOURCE_DIR}/CMake/config.h.in
               ${CMAKE_BINARY_DIR}/config.h)
configure_file(${CMAKE_SOURCE_DIR}/CMake/version.h.in
               ${CMAKE_BINARY_DIR}/version.h)

install(FILES ${CMAKE_BINARY_DIR}/config.h ${CMAKE_BINARY_DIR}/version.h
        DESTINATION "include/${CMAKE_PROJECT_NAME}")

#------------------------------------------------------------------------------
if(BUILD_WITH_PYTHON)
  add_subdirectory(python)
endif()

#------------------------------------------------------------------------------
# Documentation
find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(
    doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating developer API documentation with Doxygen"
    VERBATIM)
elseif(DOXYGEN_FOUND)

endif(DOXYGEN_FOUND)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Allow packaging with "make package"
if(BUILD_APT_PACKAGES)
  add_subdirectory(package)
endif()
#------------------------------------------------------------------------------
