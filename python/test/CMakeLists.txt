# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# List of tests
set(TEST tpybindings)

file(GLOB PYTESTFILES "test_*.py")

# Add tests to "pybindings", can be invoked with
# "ctest -L pybindings"
add_test(NAME ${TEST} COMMAND python3 -m pytest -v -s --exitfirst
                              --junit-xml=${TEST}.xml ${PYTESTFILES})
set_tests_properties(
  ${TEST}
  PROPERTIES
    LABELS
    "pybindings"
    # Make sure mocks are downloaded before invoking
    # pybindings test (to avoid race conditions)
    DEPENDS
    download_mocks
    ENVIRONMENT
    "PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH};DATA_DIR=${DATA_DIR};SCRIPTS_DIR=${CMAKE_SOURCE_DIR}/scripts"
)
