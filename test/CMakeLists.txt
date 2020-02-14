set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures)

set(TEST_MEASUREMENTSET CACHE STRING "measurement set used for testing")

if (TEST_MEASUREMENTSET)

  # Find legacy BeamModel library (e.g. version 4.0)
  find_package(BeamModel)

  if (BEAMMODEL_FOUND)
  # Test with legacy beam model library
  add_executable(tStationReference tStation.cc)
  target_link_libraries(
    tStationReference PRIVATE
    ${BEAMMODEL_STATION_RESPONSE_LIBRARY}
    ${CASACORE_LIBRARIES}
  )
  target_include_directories(
    tStationReference PRIVATE
    ${BEAMMODEL_INCLUDE_DIR}/StationResponse
    ${CASACORE_INCLUDE_DIR}
    "${CMAKE_SOURCE_DIR}"
    "${CMAKE_BINARY_DIR}"
  )
  add_test(station tStationReference)
  endif()

  # Test with this beam model library
  add_executable(tStation tStation.cc)
  target_link_libraries(
    tStation PRIVATE
    stationresponse
  )
  target_include_directories(
    tStation PRIVATE
    "${CMAKE_SOURCE_DIR}"
    "${CMAKE_BINARY_DIR}"
  )
  add_test(station tStation)
endif()

find_package(CFITSIO)

if (CFITSIO_FOUND)
  add_library(fitsio SHARED fitsreader.cpp fitswriter.cpp fitsiochecker.cpp)
  target_include_directories(fitsio PUBLIC ${CFITSIO_INCLUDE_DIR})
  target_include_directories(fitsio PUBLIC ${CASACORE_INCLUDE_DIR})
  target_link_libraries(fitsio PUBLIC ${CFITSIO_LIBRARY})
  target_link_libraries(fitsio PUBLIC ${CASACORE_LIBRARIES})
  target_link_libraries(fitsio PUBLIC ${CASA_FITS_LIBRARY}) # TODO: why is this needed?

  add_library(beamhelper SHARED beam-helper.cpp)
  target_include_directories(beamhelper PUBLIC "${CMAKE_SOURCE_DIR}")
  target_link_libraries(beamhelper PUBLIC stationresponse)
  target_link_libraries(beamhelper PUBLIC fitsio)

  add_executable(tStationBeamHamaker.x tStationBeamHamaker.cc)
  target_link_libraries(tStationBeamHamaker.x PUBLIC beamhelper)
  target_include_directories(tStationBeamHamaker.x PUBLIC ${CMAKE_BINARY_DIR})

  add_executable(tStationBeamOSKARDipole.x tStationBeamOSKARDipole.cc)
  target_link_libraries(tStationBeamOSKARDipole.x PUBLIC beamhelper)
  target_include_directories(tStationBeamOSKARDipole.x PUBLIC ${CMAKE_BINARY_DIR})

  add_executable(tStationBeamOSKARSphericalWave.x tStationBeamOSKARSphericalWave.cc)
  target_link_libraries(tStationBeamOSKARSphericalWave.x PUBLIC beamhelper)
  target_include_directories(tStationBeamOSKARSphericalWave.x PUBLIC ${CMAKE_BINARY_DIR})

  find_package(OpenMP)
  if (OpenMP_CXX_FOUND)
      target_link_libraries(tStationBeamHamaker.x PUBLIC OpenMP::OpenMP_CXX)
      target_link_libraries(tStationBeamOSKARDipole.x PUBLIC OpenMP::OpenMP_CXX)
      target_link_libraries(tStationBeamOSKARSphericalWave.x PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()