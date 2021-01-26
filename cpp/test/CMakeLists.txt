# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Find boost
find_package(Boost COMPONENTS unit_test_framework)

# Set cache values if not specified
set(LOFAR_HBA_MOCK_MS ${CMAKE_SOURCE_DIR}/test_data/LOFAR_HBA_MOCK.ms CACHE PATH "")
set(LOFAR_LBA_MOCK_MS ${CMAKE_SOURCE_DIR}/test_data/LOFAR_LBA_MOCK.ms CACHE PATH "")
set(OSKAR_MOCK_MS ${CMAKE_SOURCE_DIR}/test_data/OSKAR_MOCK.ms CACHE PATH "")
set(VLA_MOCK_MS ${CMAKE_SOURCE_DIR}/test_data/VLA_MOCK.ms CACHE PATH "")
set(MWA_MOCK_MS ${CMAKE_SOURCE_DIR}/test_data/MWA_MOCK.ms CACHE PATH "")
set(MWA_COEFF_PATH ${CMAKE_SOURCE_DIR}/test_data/mwa_full_embedded_element_pattern.h5 CACHE PATH "")
set(MOCK_H5PARM_PATH ${CMAKE_SOURCE_DIR}/test_data/MOCK_H5PARM.h5 CACHE PATH "")

#------------------------------------------------------------------------------
set(TEST_FILENAMES
  runtests.cc
  tlofar_hba.cc
  tlofar_lba.cc
  toskar.cc
  tmwa.cc
  tstation.cc
  tvla.cc
  tcache.cc
  th5parmaterm.cc
)

# Add boost dynamic link flag for all test files.
# https://www.boost.org/doc/libs/1_66_0/libs/test/doc/html/boost_test/usage_variants.html
# Without this flag, linking is incorrect and boost performs duplicate delete()
# calls after running all tests, in the cleanup phase.
set_source_files_properties(
  ${TEST_FILENAMES} PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK"
)

add_test(
  NAME download_lofar_hba_ms
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_lofar_hba_ms.sh
)
add_test(
  NAME download_lofar_lba_ms
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_lofar_lba_ms.sh
)
add_test(
  NAME download_oskar_ms
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_oskar_ms.sh
)
add_test(
  NAME download_vla_ms
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_vla_ms.sh
)
add_test(
  NAME download_mwa_ms
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_ms.sh
)
add_test(
  NAME download_mwa_coeff
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_coeff.sh
)
set_tests_properties(
  download_lofar_hba_ms download_lofar_lba_ms download_oskar_ms download_vla_ms download_mwa_ms download_mwa_coeff
  PROPERTIES FIXTURES_SETUP mocks
)
add_test(
  NAME download_mock_h5parm
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mock_h5parm.sh
)
set_tests_properties(
  download_lofar_hba_ms download_lofar_lba_ms download_oskar_ms download_vla_ms download_mwa_ms download_mwa_coeff download_mock_h5parm
  PROPERTIES FIXTURES_SETUP mocks
)

add_executable(
  unittests EXCLUDE_FROM_ALL
  ${TEST_FILENAMES}
)
target_link_libraries(unittests everybeam ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} schaapcommon)
# Required to get the config.h header
target_include_directories(unittests PRIVATE "${CMAKE_BINARY_DIR}")

# Add tests for ctest
add_test(buildunittests ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target unittests)
set_tests_properties(buildunittests PROPERTIES FIXTURES_SETUP unittests)
add_test(
  NAME unittests
  COMMAND unittests -f JUNIT -k unittests.xml --catch_system_error=yes
)
set_tests_properties(
  unittests PROPERTIES LABELS unit
  FIXTURES_REQUIRED "unittests;mocks"
)
