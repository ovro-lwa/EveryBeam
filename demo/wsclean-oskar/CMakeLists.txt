# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

set(MS "testdata.ms")

# Generate the test MeasurementSet
add_custom_command(
  OUTPUT testdata.ms
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND
    ${CMAKE_COMMAND} -E env PATH="${CMAKE_SOURCE_DIR}/scripts/misc:$ENV{PATH}"
    PYTHONPATH="${CMAKE_CURRENT_SOURCE_DIR}/../utils:$ENV{PYTHONPATH}"
    DATA_DIR="${DATA_DIR}" python3 -B
    "${CMAKE_CURRENT_SOURCE_DIR}/run_oskar_simulation.py")

# Generate a template image using the pointing in the test MeasurementSet
add_custom_command(
  OUTPUT template-image.fits
  COMMAND wsclean -name template -interval 0 1 -no-dirty -scale 2amin -size 300
          300 ${MS}
  DEPENDS testdata.ms)

set(MODEL_IMAGES wsclean-I-model.fits wsclean-Q-model.fits wsclean-U-model.fits
                 wsclean-V-model.fits)

# Use template image to generate model images with pointsources
add_custom_command(
  OUTPUT ${MODEL_IMAGES}
  COMMAND python3 -B "${CMAKE_CURRENT_SOURCE_DIR}/make_model_images.py"
  DEPENDS template-image.fits make_model_images.py)

# Main target to run the wsclean-oskar tests
# Uses the test MeasurementSet and model images as input
add_custom_target(
  run-wsclean-oskar
  COMMAND
    ${CMAKE_COMMAND} -E env
    EXTRA_PATH="${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_SOURCE_DIR}/scripts/coeff_scripts"
    DATA_DIR="${DATA_DIR}" MS=${MS}
    "${CMAKE_CURRENT_SOURCE_DIR}/run_wsclean_oskar.sh"
  DEPENDS ${MODEL_IMAGES} testdata.ms run_wsclean_oskar.sh)
