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
set(CMAKE_INSTALL_DATA_DIR "${CMAKE_INSTALL_PREFIX}/share/${projectname}")
message("Storing data in: " ${CMAKE_INSTALL_DATA_DIR})
configure_file(${CMAKE_SOURCE_DIR}/CMake/config.h.in ${CMAKE_BINARY_DIR}/config.h)

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

#------------------------------------------------------------------------------
# Set CMake and compiler options
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()

if (POLICY CMP0074)
cmake_policy(SET CMP0074 NEW)
endif()

# Set compile options
add_compile_options(-std=c++11 "${OpenMP_CXX_FLAGS}" -Wall -DNDEBUG -Wl,--no-undefined)

#------------------------------------------------------------------------------
# Find/load top level dependencies
add_subdirectory(external)

#------------------------------------------------------------------------------
# Add source
add_subdirectory(cpp)

#------------------------------------------------------------------------------
# Add tests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
   add_subdirectory(cpp/test)
   # TODO: compiling the demos should probably be a different cmake project
   # in which we use find_package(EveryBeam)
   add_subdirectory(demo)
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
