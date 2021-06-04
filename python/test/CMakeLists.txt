# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# List of tests
set(UNIT_TESTS tPyLoad)

# Add tests to "pybindings", can be invoked with
# "ctest -L pybindings"
foreach(TEST ${UNIT_TESTS})
  add_test(${TEST} ${CMAKE_CURRENT_SOURCE_DIR}/${TEST}.sh)
  set_tests_properties(
    ${TEST}
    PROPERTIES
      LABELS
      "pybindings"
      ENVIRONMENT
      "LIB_DIR=${CMAKE_BINARY_DIR}/python;DATA_DIR=${CMAKE_SOURCE_DIR}/test_data;SCRIPTS_DIR=${CMAKE_SOURCE_DIR}/scripts"
  )
endforeach()
