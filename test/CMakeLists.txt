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
