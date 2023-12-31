# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# directory for config.h
include_directories(${CMAKE_BINARY_DIR})

# build libhamaker.so
add_library(hamaker SHARED hamakerelementresponse.cc hamakercoeff.cc)
target_link_libraries(hamaker PUBLIC ${HDF5_CXX_LIBRARIES} everybeam-core)
string(TOLOWER ${CMAKE_PROJECT_NAME} projectname)
set_target_properties(hamaker PROPERTIES LIBRARY_OUTPUT_NAME
                                         "${projectname}-hamaker")

# install libhamaker.so
install(
  TARGETS hamaker
  EXPORT EveryBeamTargets
  DESTINATION lib)

# install coefficients
message("install hamaker in: " ${EVERYBEAM_FULL_DATADIR})
install(FILES "${CMAKE_SOURCE_DIR}/coeffs/HamakerHBACoeff.h5"
        DESTINATION ${EVERYBEAM_DATADIR})
install(FILES "${CMAKE_SOURCE_DIR}/coeffs/HamakerLBACoeff.h5"
        DESTINATION ${EVERYBEAM_DATADIR})
