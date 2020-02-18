set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures)


add_executable(tStationReference tStation.cc)
target_link_libraries(tStationReference "/home/vdtol/src/LOFARBeam/build/install/lib/libstationresponse.so" ${CASACORE_LIBRARIES})
target_include_directories(tStationReference PRIVATE
  "/home/vdtol/src/LOFARBeam/build/install/include/StationResponse"
  ${CASACORE_INCLUDE_DIR}
  "${CMAKE_SOURCE_DIR}"
)

add_executable(tStation tStation.cc)
target_link_libraries(tStation stationresponse)
target_include_directories(tStation PRIVATE "${CMAKE_SOURCE_DIR}")



add_test(station tStation)
