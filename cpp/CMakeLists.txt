# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

project(everybeam_cpp)

#------------------------------------------------------------------------------
add_subdirectory(aterms)
add_subdirectory(common)
add_subdirectory(coords)
add_subdirectory(pointresponse)
add_subdirectory(griddedresponse)
add_subdirectory(hamaker)
add_subdirectory(oskar)
add_subdirectory(telescope)
add_subdirectory(circularsymmetric)
add_subdirectory(lobes)

#------------------------------------------------------------------------------
add_library(
  everybeam SHARED
  antenna.cc
  elementresponse.cc
  beamformer.cc
  beamformeridenticalantennas.cc
  beamformerlofar.cc
  beamformerlofarhba.cc
  beamformerlofarlba.cc
  element.cc
  elementhamaker.cc
  load.cc
  # aterms
  aterms/atermconfig.cc
  aterms/atermbase.cc
  aterms/atermresampler.cc
  aterms/everybeamaterm.cc
  aterms/fitsatermbase.cc
  aterms/fitsaterm.cc
  aterms/dldmaterm.cc
  aterms/pafbeamterm.cc
  aterms/h5parmaterm.cc
  aterms/cache.cc
  # common
  common/sphericalharmonics.cc
  common/fftresampler.cc
  coords/itrfconverter.cc
  coords/itrfdirection.cc
  msreadutils.cc
  station.cc
  phasedarrayresponse.cc
  telescope/lofar.cc
  telescope/dish.cc
  telescope/mwa.cc
  telescope/phasedarray.cc
  telescope/oskar.cc
  griddedresponse/griddedresponse.cc
  griddedresponse/dishgrid.cc
  griddedresponse/mwagrid.cc
  griddedresponse/aartfaacgrid.cc
  pointresponse/phasedarraypoint.cc
  pointresponse/dishpoint.cc
  pointresponse/mwapoint.cc
  circularsymmetric/atcacoefficients.cc
  circularsymmetric/gmrtcoefficients.cc
  circularsymmetric/vlacoefficients.cc
  circularsymmetric/voltagepattern.cc
  # Phased array telescopes (SKA, LOFAR)
  griddedresponse/phasedarraygrid.cc
  # MWA(beam) related
  mwabeam/tilebeam2016.cc
  mwabeam/beam2016implementation.cc)

# Make sure that when other targets within this project link against the everybeam target,
# they can find the include files.
target_include_directories(
  everybeam PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
# Required to get config.h header
target_include_directories(everybeam PRIVATE "${CMAKE_BINARY_DIR}")

target_include_directories(everybeam PUBLIC ${CASACORE_INCLUDE_DIR})
target_link_libraries(everybeam PUBLIC hamaker oskar lobes)
target_link_libraries(
  everybeam PUBLIC ${CASACORE_LIBRARIES} ${HDF5_CXX_LIBRARIES}
                   ${FFTW3F_LIBRARY} ${CFITSIO_LIBRARY})
target_link_libraries(everybeam PRIVATE schaapcommon)

install(
  TARGETS everybeam
  EXPORT EveryBeamTargets
  DESTINATION lib)

install(
  FILES antenna.h
        beamformer.h
        beamformeridenticalantennas.h
        beammode.h
        beamnormalisationmode.h
        correctionmode.h
        element.h
        elementhamaker.h
        elementresponse.h
        fieldresponse.h
        msreadutils.h
        station.h
        load.h
        options.h
        phasedarrayresponse.h
  DESTINATION "include/${CMAKE_PROJECT_NAME}")

# TODO: not sure whether this is needed at all
install(
  EXPORT EveryBeamTargets
  FILE EveryBeamTargets.cmake
  NAMESPACE EveryBeam::
  DESTINATION "lib/${projectname}")

# Set-up cmake configuration files
configure_file("${CMAKE_SOURCE_DIR}/CMake/EveryBeamConfig.cmake.in"
               "${CMAKE_BINARY_DIR}/CMakeFiles/EveryBeamConfig.cmake" @ONLY)
configure_file(
  "${CMAKE_SOURCE_DIR}/CMake/EveryBeamConfigVersion.cmake.in"
  "${CMAKE_BINARY_DIR}/CMakeFiles/EveryBeamConfigVersion.cmake" @ONLY)

install(FILES "${CMAKE_BINARY_DIR}/CMakeFiles/EveryBeamConfig.cmake"
              "${CMAKE_BINARY_DIR}/CMakeFiles/EveryBeamConfigVersion.cmake"
        DESTINATION "lib/${projectname}")
