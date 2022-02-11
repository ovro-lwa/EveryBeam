# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

#------------------------------------------------------------------------------
# directory for config.h
include_directories(${CMAKE_BINARY_DIR})

#------------------------------------------------------------------------------
add_library(lobes SHARED lobeselementresponse.cc
                         ${CMAKE_SOURCE_DIR}/cpp/common/sphericalharmonics.cc)
target_include_directories(lobes PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(
  lobes PUBLIC ${HDF5_CXX_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY}
               ${Boost_SYSTEM_LIBRARY})
string(TOLOWER ${CMAKE_PROJECT_NAME} projectname)
set_target_properties(lobes PROPERTIES LIBRARY_OUTPUT_NAME
                                       "${projectname}-lobes")

#------------------------------------------------------------------------------
if(DOWNLOAD_LOBES)
  # Download the LOBES coefficient files if DOWNLOAD_LOBES=On, leads to a relatively small amount of overhead
  # in make command, even if the coefficient files have been dowloaded
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/coeffs/lobes)
  add_custom_target(
    download_lobes_coefficients
    COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_lobes_coeffs.sh
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/coeffs/lobes)

  #Execute downloading of LOBEs coefficients here
  add_dependencies(lobes download_lobes_coefficients)
endif()

#------------------------------------------------------------------------------
# Install
install(
  TARGETS lobes
  EXPORT EveryBeamTargets
  DESTINATION lib)

# install coefficients
message("install lobes coefficients in: " ${EVERYBEAM_FULL_DATADIR}/lobes)
install(
  DIRECTORY "${CMAKE_BINARY_DIR}/coeffs/lobes"
  DESTINATION ${EVERYBEAM_DATADIR}
  FILES_MATCHING
  PATTERN "LOBES_*")

#------------------------------------------------------------------------------
# TODO: can't we remove all this?
set(MISC_DIR ${CMAKE_SOURCE_DIR}/scripts/misc)

configure_file(${MISC_DIR}/F4far_new.py
               ${CMAKE_CURRENT_BINARY_DIR}/F4far_new.py COPYONLY)
configure_file(${MISC_DIR}/calc_modes_1deg_noback.py
               ${CMAKE_CURRENT_BINARY_DIR}/calc_modes_1deg_noback.py COPYONLY)
configure_file(${MISC_DIR}/hamaker_vs_lobes.py
               ${CMAKE_CURRENT_BINARY_DIR}/hamaker_vs_lobes.py COPYONLY)
configure_file(${MISC_DIR}/test_beam_model.py
               ${CMAKE_CURRENT_BINARY_DIR}/test_beam_model.py COPYONLY)
