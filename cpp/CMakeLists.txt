#------------------------------------------------------------------------------
add_subdirectory(common)
add_subdirectory(coords)
add_subdirectory(gridded_response)
add_subdirectory(hamaker)
add_subdirectory(lobes)
add_subdirectory(oskar)
add_subdirectory(telescope)

#------------------------------------------------------------------------------
add_library(everybeam SHARED
  Antenna.cc
  ElementResponse.cc
  BeamFormer.cc
  Element.cc
  load.cc
  coords/ITRFConverter.cc
  coords/ITRFDirection.cc
  LofarMetaDataUtil.cc
  Station.cc
  telescope/lofar.cc
  gridded_response/lofargrid.cc
)

target_include_directories(everybeam PUBLIC ${CASACORE_INCLUDE_DIR})
target_link_libraries(everybeam PUBLIC ${CASACORE_LIBRARIES})
target_link_libraries(everybeam PUBLIC hamaker lobes oskar)

install (
  TARGETS everybeam
  EXPORT EveryBeamTargets
  DESTINATION lib)

install (FILES
  Antenna.h
  BeamFormer.h
  Element.h
  ElementResponse.h
  LofarMetaDataUtil.h
  Station.h
  # Related to new API:
  load.h
  options.h
DESTINATION "include/${CMAKE_PROJECT_NAME}")

install(
  EXPORT EveryBeamTargets
  FILE EveryBeamConfig.cmake
  NAMESPACE EveryBeam::
  DESTINATION "lib/${projectname}"
)

