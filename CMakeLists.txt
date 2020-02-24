cmake_minimum_required(VERSION 3.0.0)

project(lofarbeam)

# Configure directory for data files
set(CMAKE_INSTALL_DATA_DIR "${CMAKE_INSTALL_PREFIX}/data")
message("Storing data in: " ${CMAKE_INSTALL_DATA_DIR})
configure_file(${CMAKE_SOURCE_DIR}/CMake/config.h.in ${CMAKE_BINARY_DIR}/config.h)

add_subdirectory(external)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()

# Add tests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(test)
endif()

if (POLICY CMP0074)
cmake_policy(SET CMP0074 NEW)
endif()

find_package(HDF5 COMPONENTS C CXX REQUIRED)
add_definitions(${HDF5_DEFINITIONS})
include_directories(${HDF5_INCLUDE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)
include_directories(${CASACORE_INCLUDE_DIR})

# Add element models
add_subdirectory(hamaker)
add_subdirectory(lobes)
add_subdirectory(oskar)

add_compile_options(-std=c++11 -Wall -DNDEBUG -Wl,--no-undefined)

add_library(stationresponse SHARED
  Antenna.cc
  ElementResponse.cc
  BeamFormer.cc
  Element.cc
  ITRFConverter.cc
  ITRFDirection.cc
  LofarMetaDataUtil.cc
  MathUtil.cc
  Station.cc
  Types.cc
)

set_target_properties(stationresponse PROPERTIES
    CXX_STANDARD 17
    VERSION 3
)

target_include_directories(stationresponse PUBLIC ${CASACORE_INCLUDE_DIR})

target_link_libraries(stationresponse PUBLIC ${CASACORE_LIBRARIES})

target_link_libraries(stationresponse PUBLIC hamaker lobes oskar)

install (TARGETS stationresponse DESTINATION lib)

install (FILES
  Antenna.h
  BeamFormer.h
  Constants.h
  Element.h
  ElementResponse.h
  LofarMetaDataUtil.h
  MutablePtr.h
  Station.h
  Types.h
  ITRFConverter.h
  ITRFDirection.h
DESTINATION include/StationResponse)

find_package(Doxygen)

if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating developer API documentation with Doxygen" VERBATIM)
elseif(DOXYGEN_FOUND)

endif(DOXYGEN_FOUND)
