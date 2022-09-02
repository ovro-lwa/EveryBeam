# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Ensure that the Boost unit test framework is available.
find_package(
  Boost
  COMPONENTS unit_test_framework
  REQUIRED)

# Set cache values if not specified
set(MOCK_MS
    LOFAR_HBA_MOCK
    LOFAR_LBA_MOCK
    AARTFAAC_LBA_MOCK
    AARTFAAC_6_LBA_MOCK
    OSKAR_MOCK
    SCREEN_FITTING
    SKA_MID_MOCK
    VLA_MOCK
    MWA_MOCK)
foreach(mock ${MOCK_MS})
  set(${mock}_MS
      ${DATA_DIR}/${mock}.ms
      CACHE PATH "")
endforeach()

set(MWA_COEFF_PATH
    ${DATA_DIR}/mwa_full_embedded_element_pattern.h5
    CACHE PATH "")
set(MOCK_H5PARM_PATH
    ${DATA_DIR}/MOCK_H5PARM.h5
    CACHE PATH "")
set(H5_SOLUTIONS_PATH
    ${DATA_DIR}/solutions.h5
    CACHE PATH "")
set(KL_SCREEN_PATH
    ${DATA_DIR}/kl_screen_0.fits
    CACHE PATH "")
set(LOBES_COEFF_PATH
    ${CMAKE_BINARY_DIR}/coeffs/lobes
    CACHE PATH "")

#------------------------------------------------------------------------------
set(TEST_FILENAMES
    runtests.cc
    tatermconfig.cc
    tgriddedresponse.cc
    tlofar_hba.cc
    tlofar_lba.cc
    tlofar_aartfaac.cc
    toskar.cc
    tmwa.cc
    tphasedarraypoint.cc
    tskamid.cc
    tstation.cc
    tvla.cc
    tcache.cc
    tklfittingaterm.cc
    tfourierfittingaterm.cc
    tfourierfitter.cc
    tklfitter.cc
    th5parmaterm.cc)

# Add boost dynamic link flag for all test files.
# https://www.boost.org/doc/libs/1_66_0/libs/test/doc/html/boost_test/usage_variants.html
# Without this flag, linking is incorrect and boost performs duplicate delete()
# calls after running all tests, in the cleanup phase.
set_source_files_properties(${TEST_FILENAMES} PROPERTIES COMPILE_DEFINITIONS
                                                         "BOOST_TEST_DYN_LINK")

ExternalProject_Add(
  download_mock_aartfaac_6_lba
  EXCLUDE_FROM_ALL TRUE
  URL https://support.astron.nl/software/ci_data/EveryBeam/AARTFAAC_6_LBA_MOCK.ms.tar.gz
  URL_HASH
    SHA256=8d7f93ce089b3a29c4d0ee8b367d27a821175d0626793b16be1643ee9bcdc6a5
  SOURCE_DIR ${DATA_DIR}/AARTFAAC_6_LBA_MOCK.ms
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")

add_custom_target(
  download_mocks
  WORKING_DIRECTORY ${DATA_DIR}
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_lofar_hba_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_lofar_lba_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_aartfaac_lba_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_h5_solutions.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_oskar_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_ska_mid_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_vla_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_ms.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mwa_coeff.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_mock_h5parm.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_kl_fits_screen.sh
  COMMAND ${CMAKE_SOURCE_DIR}/scripts/download_screen_fitting_ms.sh
  DEPENDS download_mock_aartfaac_6_lba)
add_test(NAME download_mocks
         COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target
                 download_mocks)
set_tests_properties(download_mocks PROPERTIES FIXTURES_SETUP mocks)

add_executable(unittests ${TEST_FILENAMES})
target_link_libraries(unittests everybeam ${Boost_LIBRARIES} schaapcommon)

target_include_directories(
  unittests
  PRIVATE "${CMAKE_BINARY_DIR}" # Required to get the config.h header
          "${CMAKE_SOURCE_DIR}/external/xtensor/include/"
          "${CMAKE_SOURCE_DIR}/external/xtensor-fftw/include/"
          "${CMAKE_SOURCE_DIR}/external/xtl/include/")

set_property(TARGET unittests PROPERTY CXX_STANDARD 17)

# Add test for automatically (re)building unittests if needed.
add_test(NAME buildunittests COMMAND ${CMAKE_COMMAND} --build
                                     ${CMAKE_BINARY_DIR} --target unittests)
set_tests_properties(buildunittests PROPERTIES FIXTURES_SETUP unittests)

add_test(NAME unittests COMMAND unittests -f JUNIT -k unittests.xml
                                --catch_system_error=yes)
set_tests_properties(unittests PROPERTIES LABELS unit FIXTURES_REQUIRED
                                          "unittests;mocks")
