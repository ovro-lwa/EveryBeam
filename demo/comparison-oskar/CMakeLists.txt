# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

#------------------------------------------------------------------------------
# CMake file for compiling a comparison between OSKAR and EveryBeam
add_executable(make_element_response_image make_element_response_image.cpp)
target_link_libraries(make_element_response_image oskar)

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/telescope.tm/layout.txt"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/telescope.tm")

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/telescope.tm/position.txt"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/telescope.tm")

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/telescope.tm/station000/layout.txt"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/station000/telescope.tm")

#------------------------------------------------------------------------------
# comparison-oskar knits together the cpp code and the python scripts
add_custom_target(
  comparison-oskar-basefunctions
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env
    PATH="${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_SOURCE_DIR}/scripts/coeff_scripts:$ENV{PATH}"
    PYTHONPATH="${CMAKE_CURRENT_SOURCE_DIR}/../utils:$ENV{PYTHONPATH}" python3
    -B "${CMAKE_CURRENT_SOURCE_DIR}/generate_basefunction_plots.py"
  DEPENDS make_element_response_image)

add_custom_target(
  comparison-oskar-station-response
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env
    PATH="${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_SOURCE_DIR}/scripts/coeff_scripts:${CMAKE_SOURCE_DIR}/scripts/misc:$ENV{PATH}"
    PYTHONPATH="${CMAKE_CURRENT_SOURCE_DIR}/../utils:$ENV{PYTHONPATH}"
    DATA_DIR=${DATA_DIR} python3 -B
    "${CMAKE_CURRENT_SOURCE_DIR}/compare_stationresponse.py")
