# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# EveryBeam uses a hard-coded, compile-time generated, path to the location
# of the data files it uses. In order to avoid that the data files in the
# package are installed in a different directory than where they are searched
# for, we need to follow the advice given in the merge request
# "Set default value of CPACK_PACKAGING_INSTALL_PREFIX properly"
# (https://gitlab.kitware.com/cmake/cmake/-/merge_requests/827).
# Set `CPACK_PACKAGING_INSTALL_PREFIX` (the install prefix used by `cpack`),
# to the value of `CMAKE_INSTALL_PREFIX` (the install prefix used by `cmake`),
# unless the user specified `CPACK_PACKAGING_INSTALL_PREFIX` explicitly.
if(NOT CPACK_PACKAGING_INSTALL_PREFIX)
  set(CPACK_PACKAGING_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
endif()

# By default, Debian packages are installed in `/usr`. Print a warning when
# `CMAKE_INSTALL_PREFIX` is set to a different directory, because that might
# break compatibility with other Debian packages.
if(NOT "${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr")
  message(
    WARNING
      "[CPack] - CMAKE_INSTALL_PREFIX is set to `${CMAKE_INSTALL_PREFIX}`. "
      "For packaging, it is strongly recommended to set it to `/usr` "
      "to avoid compatibility issues with other Debian packages.")
endif()

# Different settings for `CPACK_PACKAGING_INSTALL_PREFIX` and
# `CMAKE_INSTALL_PREFIX` may break the package being built.
if(NOT "${CMAKE_INSTALL_PREFIX}" STREQUAL "${CPACK_PACKAGING_INSTALL_PREFIX}")
  message(
    WARNING
      "[CPack] - CMAKE_INSTALL_PREFIX is set to `${CMAKE_INSTALL_PREFIX}`, and "
      "CPACK_PACKAGING_INSTALL_PREFIX is set to `${CPACK_PACKAGING_INSTALL_PREFIX}`. "
      "For packaging, it is strongly recommend to use the same setting for both, "
      "preferably `/usr` to avoid compatibility issues with other Debian packages."
  )
endif()

# Using absolute paths when packaging may break relocation. Issue a warning
# when `EVERYBEAM_DATADIR` contains an absolute path. Also, let CPack itself
# generate a warning if it encounters (other) absolute install paths.
if(IS_ABSOLUTE ${EVERYBEAM_DATADIR})
  message(WARNING "[CPack] - EVERYBEAM_DATADIR contains an absolute path. "
                  "This may break relocation during packaging.")
endif()
set(CPACK_WARN_ON_ABSOLUTE_INSTALL_DESTINATION ON)

# We use git to determine the version string from the latest tag.
find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND bash -c "${GIT_EXECUTABLE} describe --tags --dirty"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE VERSION_STRING
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  # Strip leading non-numeric characters, like "v"
  string(REGEX REPLACE "^[^[0-9]]*" "" VERSION_STRING ${VERSION_STRING})
endif()
if(NOT VERSION_STRING)
  message(
    WARNING "[CPack] Could not find a git tag to a create version string, "
            "trying CMake variable PROJECT_VERSION instead.")
  if(PROJECT_VERSION)
    set(VERSION_STRING "${PROJECT_VERSION}")
  else()
    message(FATAL_ERROR "[CPack] Failed to create version string. Bailing out!")
  endif()
endif()

set(CPACK_PACKAGE_NAME "everybeam")
set(CPACK_PACKAGE_VENDOR "ASTRON")
set(CPACK_PACKAGE_VERSION "${VERSION_STRING}")

set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "The EveryBeam library provides antenna response patterns for several instruments, such as LOFAR (and LOBES), SKA (OSKAR), MWA, JVLA, etc."
)
set(CPACK_PACKAGE_HOMEPAGE_URL "https://git.astron.nl/RD/EveryBeam")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")

set(CPACK_GENERATOR "DEB")
set(CPACK_SOURCE_GENERATOR "TGZ")

set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
set(CPACK_DEBIAN_PACKAGE_CONFLICTS)
set(CPACK_DEBIAN_PACKAGE_DEPENDS)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "deb-packages@astron.nl")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

include(CPack)

message(STATUS "Package name: ${CPACK_PACKAGE_NAME}")
message(STATUS "Package version: ${CPACK_PACKAGE_VERSION}")
