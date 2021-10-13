# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

project(everybeam_pybind11)

find_package(PythonLibs 3 REQUIRED)
find_package(PythonInterp REQUIRED)
message(STATUS "Using python version ${PYTHON_VERSION_STRING}")

# Create the binding library
# pyeverybeam is a temporary alias for everybeam target
pybind11_add_module(
  pyeverybeam wrappers/pyutils.cc wrappers/pyload.cc wrappers/pytelescope.cc
  wrappers/pylobes.cc wrappers/wrapper.cc)

target_include_directories(pyeverybeam
                           PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/cpp>")
target_link_libraries(pyeverybeam PUBLIC everybeam)
set_target_properties(pyeverybeam PROPERTIES OUTPUT_NAME everybeam)

# If the PYTHON_LIBRARY_DIR is not specified, install in lib dir
if(NOT DEFINED ${PYTHON_LIBRARY_DIR})
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c
            "import site; print(site.getsitepackages()[0])"
    OUTPUT_VARIABLE PYTHON_DIST_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(PYTHON_DIST_PATH MATCHES
     "\\/(python${PYTHON_VERSION_MAJOR}\\.${PYTHON_VERSION_MINOR}\\/.*)")
    set(PYTHON_LIBRARY_DIR ${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_MATCH_1})
  else()
    message(
      FATAL_ERROR "Failed to parse PYTHON_DIST_PATH='${PYTHON_DIST_PATH}'")
  endif()
endif()

# Install pyeverybeam in site-packages directory
install(
  TARGETS pyeverybeam
  COMPONENT python
  LIBRARY DESTINATION ${PYTHON_LIBRARY_DIR})
