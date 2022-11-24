# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

#------------------------------------------------------------------------------
# directory for config.h
include_directories(${CMAKE_BINARY_DIR})

#------------------------------------------------------------------------------
# build libeverybeam-oskar.so
add_library(oskar SHARED oskarelementresponse.cc oskardatafile.cc
                         oskardataset.cc)

string(TOLOWER ${CMAKE_PROJECT_NAME} projectname)
set_target_properties(oskar PROPERTIES LIBRARY_OUTPUT_NAME
                                       "${projectname}-oskar")

# Make sure that when other targets within this project link against the oskar target,
# they can find the include files.
target_include_directories(
  oskar PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

#------------------------------------------------------------------------------
# Link against HDF5 and OpenMP
target_link_libraries(oskar PUBLIC ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES}
                                   ${OpenMP_CXX_FLAGS} everybeam-core)
target_link_libraries(oskar PRIVATE oskar_beam_utils)

#------------------------------------------------------------------------------
# install libeverybeam-oskar.so
install(
  TARGETS oskar
  EXPORT EveryBeamTargets
  DESTINATION lib)

#install oskar coefficients
install(FILES "${CMAKE_SOURCE_DIR}/coeffs/oskar.h5"
        DESTINATION ${EVERYBEAM_DATADIR})
