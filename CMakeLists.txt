#------------------------------------------------------------------------------
# Top level CMakeLists.txt file for EveryBeam
cmake_minimum_required(VERSION 3.0.0)

#------------------------------------------------------------------------------
# Set project name and version number
project(EveryBeam VERSION 5.0.0)

string(TOLOWER ${CMAKE_PROJECT_NAME} projectname )

# Set the path to CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

# Configure directory for data files
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -DNDEBUG")
set(CMAKE_INSTALL_DATA_DIR "${CMAKE_INSTALL_PREFIX}/share/${projectname}")
message("Storing data in: " ${CMAKE_INSTALL_DATA_DIR})

# Find and include git submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --checkout
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

# Include aocommon/eigen3/pybind11 headers
include_directories("${CMAKE_SOURCE_DIR}/external/aocommon/include/")
include_directories("${CMAKE_SOURCE_DIR}/external/eigen/")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/pybind11")

# Find and include HDF5
find_package(HDF5 COMPONENTS C CXX REQUIRED)
add_definitions(${HDF5_DEFINITIONS})
include_directories(${HDF5_INCLUDE_DIR})

# Find and include Casacore
set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)
include_directories(${CASACORE_INCLUDE_DIR})

# Find and include OpenMP
find_package(OpenMP REQUIRED)

# Find and include Boost headers (boost::math required for MWA beam)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Find and include FFTW3 float libraries
find_package(FFTW3 COMPONENTS single REQUIRED)

#------------------------------------------------------------------------------
# Set CMake and compiler options
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()

if (POLICY CMP0074)
cmake_policy(SET CMP0074 NEW)
endif()

# Set compile options
add_compile_options(-std=c++11 "${OpenMP_CXX_FLAGS}" -Wall)
string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")

if (NOT CMAKE_BUILD_TYPE MATCHES Debug)
    add_compile_options(-DNDEBUG)
endif()

# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

#------------------------------------------------------------------------------
# Add source
add_subdirectory(cpp)

#------------------------------------------------------------------------------
# Add tests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  # TODO: compiling the demos should probably be a different cmake project
  # in which we use find_package(EveryBeam)
  add_subdirectory(demo)

  if(Boost_FOUND)
    add_subdirectory(cpp/test)
  else()
    message("Boost NOT found. This is not an error but indicates that the tests won't be built.")
  endif()
endif()

#------------------------------------------------------------------------------
# Documentation
find_package(Doxygen)

if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
       ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
       COMMENT "Generating developer API documentation with Doxygen" VERBATIM)
elseif(DOXYGEN_FOUND)

endif(DOXYGEN_FOUND)
