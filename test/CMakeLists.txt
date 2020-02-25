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
  add_library(beamhelper SHARED beam-helper.cpp)
  target_include_directories(beamhelper PUBLIC "${CMAKE_SOURCE_DIR}")
  target_include_directories(beamhelper PUBLIC ${CFITSIO_INCLUDE_DIR})
  target_link_libraries(beamhelper PUBLIC stationresponse)
  target_link_libraries(beamhelper PUBLIC ${CFITSIO_LIBRARY})

  add_executable(tStationBeamHamaker tStationBeamHamaker.cc)
  target_link_libraries(tStationBeamHamaker PUBLIC beamhelper)
  target_include_directories(tStationBeamHamaker PUBLIC ${CMAKE_BINARY_DIR})

  add_executable(tStationBeamOSKARDipole tStationBeamOSKARDipole.cc)
  target_link_libraries(tStationBeamOSKARDipole PUBLIC beamhelper)
  target_include_directories(tStationBeamOSKARDipole PUBLIC ${CMAKE_BINARY_DIR})

  add_executable(tStationBeamOSKARSphericalWave tStationBeamOSKARSphericalWave.cc)
  target_link_libraries(tStationBeamOSKARSphericalWave PUBLIC beamhelper)
  target_include_directories(tStationBeamOSKARSphericalWave PUBLIC ${CMAKE_BINARY_DIR})

  add_executable(tElementBeamHamaker tElementBeamHamaker.cc)
  target_link_libraries(tElementBeamHamaker PUBLIC beamhelper)
  target_include_directories(tElementBeamHamaker PUBLIC ${CMAKE_BINARY_DIR})

  add_executable(tElementBeamOSKARDipole tElementBeamOSKARDipole.cc)
  target_link_libraries(tElementBeamOSKARDipole PUBLIC beamhelper)
  target_include_directories(tElementBeamOSKARDipole PUBLIC ${CMAKE_BINARY_DIR})

  add_executable(tElementBeamOSKARSphericalWave tElementBeamOSKARSphericalWave.cc)
  target_link_libraries(tElementBeamOSKARSphericalWave PUBLIC beamhelper)
  target_include_directories(tElementBeamOSKARSphericalWave PUBLIC ${CMAKE_BINARY_DIR})

  find_package(OpenMP)
  if (OpenMP_CXX_FOUND)
      target_link_libraries(tStationBeamHamaker PUBLIC OpenMP::OpenMP_CXX)
      target_link_libraries(tStationBeamOSKARDipole PUBLIC OpenMP::OpenMP_CXX)
      target_link_libraries(tStationBeamOSKARSphericalWave PUBLIC OpenMP::OpenMP_CXX)
      target_link_libraries(tElementBeamHamaker PUBLIC OpenMP::OpenMP_CXX)
      target_link_libraries(tElementBeamOSKARDipole PUBLIC OpenMP::OpenMP_CXX)
      target_link_libraries(tElementBeamOSKARSphericalWave PUBLIC OpenMP::OpenMP_CXX)
  endif()

  add_test(stationbeam-hamaker tStationBeamHamaker)
  add_test(stationbeam-oskardipole tStationBeamOSKARDipole)
  add_test(stationbeam-oskarsphericalwave tStationBeamOSKARSphericalWave)
  add_test(elementbeam-hamaker tElementBeamHamaker)
  add_test(elementbeam-oskardipole tElementBeamOSKARDipole)
  add_test(elementbeam-oskarsphericalwave tElementBeamOSKARSphericalWave)

endif()
