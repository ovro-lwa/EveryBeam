# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

#------------------------------------------------------------------------------
# Top level CMakeLists.txt file for EveryBeam
cmake_minimum_required(VERSION 3.8)

#------------------------------------------------------------------------------
# Set version name and project number
set(EVERYBEAM_VERSION 0.5.1)
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
option(DOWNLOAD_LOBES "Download and install LOBEs coefficient files" OFF)

string(TOLOWER ${CMAKE_PROJECT_NAME} projectname)

# Set the path to CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)
include(ExternalProject)

# Configure directory for data files
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")
set(CMAKE_CXX_VISIBILITY_PRESET "hidden")
string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
# Note: Use type `STRING` here, instead of `PATH`. Because, if the user
# specified a _relative_ path on the command-line without specifying the type
# `PATH`, then the `set` command will treat the path as relative to the current
# working directory and convert it to an absolute path. We do *not* want this!
set(EVERYBEAM_DATADIR
    "share/${projectname}"
    CACHE STRING "EveryBeam data directory")
if(IS_ABSOLUTE ${EVERYBEAM_DATADIR})
  set(EVERYBEAM_FULL_DATADIR "${EVERYBEAM_DATADIR}")
else()
  set(EVERYBEAM_FULL_DATADIR "${CMAKE_INSTALL_PREFIX}/${EVERYBEAM_DATADIR}")
endif()
message("Storing data in: " ${EVERYBEAM_FULL_DATADIR})

# Find and include git submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  # Update submodules as needed
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --checkout
              --depth 1
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
include_directories("${CMAKE_SOURCE_DIR}/external/aocommon/include")
include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/external/eigen")
# Find python with cmake to convince pybind11 to use right python path
# This can go when we're using CMake >= 3.12:
# https://pybind11.readthedocs.io/en/stable/compiling.html#findpython-mode
find_package(PythonInterp 3)
add_subdirectory("${CMAKE_SOURCE_DIR}/external/pybind11")
include_directories(SYSTEM ${pybind11_INCLUDE_DIRS})

# Include schaapcommon
add_subdirectory("${CMAKE_SOURCE_DIR}/external/schaapcommon")
include_directories("${CMAKE_SOURCE_DIR}/external/schaapcommon/include")

# Include oskar beam utils
add_subdirectory("${CMAKE_SOURCE_DIR}/external/oskar/oskar/beam_utils")

# Find and include HDF5
find_package(
  HDF5
  COMPONENTS C CXX
  REQUIRED)
add_definitions(${HDF5_DEFINITIONS} -DH5_USE_110_API)
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})

find_package(Threads REQUIRED)

# Find and include Casacore
set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)
include_directories(SYSTEM ${CASACORE_INCLUDE_DIRS})

# Fetch XTensor and related libraries.
include(external/aocommon/CMake/FetchXTensor.cmake)

# Find and include Boost headers (boost::math required for MWA beam)
find_package(Boost REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# Find and include FFTW3 float libraries
find_package(
  FFTW3
  COMPONENTS single
  REQUIRED)
include_directories(SYSTEM ${FFTW3_INCLUDE_DIR})

# Find BLAS and LAPACK, needed for everybeam::aterms only
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Find CFITSIO, needed for everybeam::aterms only
find_package(CFITSIO REQUIRED)
include_directories(SYSTEM ${CFITSIO_INCLUDE_DIRS})

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
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
add_compile_options(
  "${OpenMP_CXX_FLAGS}"
  -Wall
  -Wnon-virtual-dtor
  -Wzero-as-null-pointer-constant
  -Wduplicated-branches
  -Wundef
  -Wvla
  -Wpointer-arith
  -Wextra
  -Wno-unused-parameter)
string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # GCC 8.x requires linking with stdc++fs for the filesystem library
  # https://gcc.gnu.org/onlinedocs/gcc-9.1.0/libstdc++/manual/manual/status.html#status.iso.2017
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    link_libraries(stdc++fs)
  elseif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0)
    message(
      FATAL_ERROR "The GCC version is too old, upgrade to GCC 8.0 or newer")
  endif()
endif()

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
# Set up a test_data directory in the build directory.
# Create symbolic links to the files/directories in the source directory.
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.14) # Required for file(CREATE_LINK)
  set(DATA_DIR ${CMAKE_BINARY_DIR}/test_data)
  file(MAKE_DIRECTORY ${DATA_DIR})
  file(GLOB SOURCE_DATA_FILES "${CMAKE_SOURCE_DIR}/test_data/*")
  foreach(SOURCE_DATA_FILE ${SOURCE_DATA_FILES})
    string(REPLACE ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} BINARY_DATA_FILE
                   ${SOURCE_DATA_FILE})
    file(CREATE_LINK ${SOURCE_DATA_FILE} ${BINARY_DATA_FILE} SYMBOLIC)
  endforeach()
else()
  # For older versions, fall back to using the source directory.
  set(DATA_DIR ${CMAKE_SOURCE_DIR}/test_data)
endif()

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

  if(BUILD_TESTING)
    add_subdirectory(python/test)
  endif()
endif()

#------------------------------------------------------------------------------
# Documentation
add_subdirectory(doc)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Allow packaging with "make package"
if(BUILD_APT_PACKAGES)
  add_subdirectory(package)
endif()
#------------------------------------------------------------------------------
