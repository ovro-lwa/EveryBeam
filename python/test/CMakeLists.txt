# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# List of tests
file(GLOB PYTESTFILES "test_*.py")

foreach(TEST ${PYTESTFILES})
  string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/test_" "py_" NAME ${TEST})
  string(REPLACE ".py" "" NAME ${NAME})

  # Add tests to "pybindings", can be invoked with
  # "ctest -L pybindings"
  add_test(NAME ${NAME} COMMAND python3 -m pytest -v -s --exitfirst
                                --junit-xml=${NAME}.xml ${TEST})
  set_tests_properties(
    ${NAME}
    PROPERTIES
      LABELS
      "pybindings"
      # Make sure mocks are downloaded before invoking
      # python test (to avoid race conditions)
      FIXTURES_REQUIRED
      mocks
      ENVIRONMENT
      "PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH};DATA_DIR=${DATA_DIR}"
  )
endforeach()
