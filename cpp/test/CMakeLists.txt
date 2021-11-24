# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Ensure that the Boost unit test framework is available.
find_package(
  Boost
  COMPONENTS unit_test_framework
  REQUIRED)

if(DOWNLOAD_LOBES)
  add_definitions(-DDOWNLOAD_LOBES)
endif()

# Set cache values if not specified
set(MOCK_MS LOFAR_HBA_MOCK LOFAR_LBA_MOCK AARTFAAC_LBA_MOCK OSKAR_MOCK VLA_MOCK
            MWA_MOCK)
foreach(mock ${MOCK_MS})
  set(${mock}_MS
      ${CMAKE_SOURCE_DIR}/test_data/${mock}.ms
      CACHE PATH "")
endforeach()

set(MWA_COEFF_PATH
    ${CMAKE_SOURCE_DIR}/test_data/mwa_full_embedded_element_pattern.h5
    CACHE PATH "")
set(MOCK_H5PARM_PATH
    ${CMAKE_SOURCE_DIR}/test_data/MOCK_H5PARM.h5
    CACHE PATH "")
set(LOBES_COEFF_PATH
    ${CMAKE_BINARY_DIR}/coeffs/lobes
    CACHE PATH "")

#------------------------------------------------------------------------------
set(TEST_FILENAMES
    runtests.cc
    tlofar_hba.cc
    tlofar_lba.cc
    tlofar_aartfaac.cc
    toskar.cc
    tmwa.cc
    tstation.cc
    tvla.cc
    tcache.cc
    th5parmaterm.cc)

# Add boost dynamic link flag for all test files.
# https://www.boost.org/doc/libs/1_66_0/libs/test/doc/html/boost_test/usage_variants.html
# Without this flag, linking is incorrect and boost performs duplicate delete()
# calls after running all tests, in the cleanup phase.
set_source_files_properties(${TEST_FILENAMES} PROPERTIES COMPILE_DEFINITIONS
                                                         "BOOST_TEST_DYN_LINK")

add_custom_target(
  download_mocks
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_lofar_hba_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_lofar_lba_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_aartfaac_lba_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_oskar_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_vla_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_coeff.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mock_h5parm.sh)
add_test(NAME download_mocks
         COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target
                 download_mocks)
set_tests_properties(download_mocks PROPERTIES FIXTURES_SETUP mocks)

add_executable(unittests ${TEST_FILENAMES})
target_link_libraries(unittests everybeam ${Boost_LIBRARIES}
                      ${CMAKE_THREAD_LIBS_INIT} schaapcommon)
# Required to get the config.h header
target_include_directories(unittests PRIVATE "${CMAKE_BINARY_DIR}")

# Add test for automatically (re)building unittests if needed.
add_test(NAME buildunittests COMMAND ${CMAKE_COMMAND} --build
                                     ${CMAKE_BINARY_DIR} --target unittests)
set_tests_properties(buildunittests PROPERTIES FIXTURES_SETUP unittests)

add_test(NAME unittests COMMAND unittests -f JUNIT -k unittests.xml
                                --catch_system_error=yes)
set_tests_properties(unittests PROPERTIES LABELS unit FIXTURES_REQUIRED
                                          "unittests;mocks")
