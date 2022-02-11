# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Ensure that we use /usr as install prefix when packaging.
set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")

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
