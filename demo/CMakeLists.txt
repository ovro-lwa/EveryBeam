# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(comparison-oskar)
add_subdirectory(wsclean-oskar)

set(TEST_MEASUREMENTSET CACHE STRING "measurement set used for testing")

if(TEST_MEASUREMENTSET)

  # Find legacy BeamModel library (e.g. version 4.0)
  find_package(BeamModel)

  if(BEAMMODEL_FOUND)
    # Test with legacy beam model library
    add_executable(tStationReference tStation.cc)
    target_link_libraries(
      tStationReference PRIVATE ${BEAMMODEL_STATION_RESPONSE_LIBRARY}
                                ${CASACORE_LIBRARIES})
    target_include_directories(
      tStationReference
      PRIVATE ${BEAMMODEL_INCLUDE_DIR}/EveryBeam ${CASACORE_INCLUDE_DIR}
              "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}")
    add_test(station tStationReference)
  endif()

  # Test with this beam model library
  add_executable(tStation tStation.cc)
  target_link_libraries(tStation PRIVATE everybeam)
  target_include_directories(tStation PRIVATE "${CMAKE_SOURCE_DIR}"
                                              "${CMAKE_BINARY_DIR}")
  add_test(station tStation)
endif()

find_package(CFITSIO)

if(CFITSIO_FOUND)
  add_library(beamhelper SHARED beam-helper.cpp)
  target_include_directories(beamhelper PUBLIC ${CMAKE_SOURCE_DIR}
                                               ${CFITSIO_INCLUDE_DIR})
  target_link_libraries(beamhelper PUBLIC everybeam ${CFITSIO_LIBRARY})

  set(CFITSIO_TESTS
      tStationBeamHamaker tStationBeamOSKARDipole
      tStationBeamOSKARSphericalWave tElementBeamHamaker
      tElementBeamOSKARDipole tElementBeamOSKARSphericalWave)
  foreach(TEST ${CFITSIO_TESTS})
    add_executable(${TEST} ${TEST}.cc)
    target_link_libraries(${TEST} PUBLIC beamhelper xtensor)
    target_include_directories(${TEST} PUBLIC ${CMAKE_BINARY_DIR})

    # These are not "tests" in the strict sense of the word, so don't add them to test stack
    # NOTE: tests should go to cpp/test
    # add_test(${TEST} ${TEST})
  endforeach()
endif()
