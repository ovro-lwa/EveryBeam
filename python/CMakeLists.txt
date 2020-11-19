# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

project(everybeam_pybind11)

find_package(PythonLibs 3 REQUIRED)
find_package(PythonInterp REQUIRED)
message(STATUS "Using python version ${PYTHON_VERSION_STRING}")

# Create the binding library
# pyeverybeam is a temporary alias for everybeam target
pybind11_add_module(pyeverybeam 
                    wrappers/pyutils.cc
                    wrappers/pyload.cc
                    wrappers/pytelescope.cc
                    wrappers/wrapper.cc
                    )


target_include_directories(pyeverybeam PUBLIC ${CASACORE_INCLUDE_DIR})
target_include_directories(pyeverybeam PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/cpp>")
target_link_libraries(pyeverybeam PUBLIC everybeam)
target_link_libraries(pyeverybeam PUBLIC ${CASACORE_LIBRARIES})
set_target_properties(pyeverybeam PROPERTIES OUTPUT_NAME everybeam)

add_subdirectory(test)

# Add the subdirectory for LOBES (convenience) functions
add_subdirectory(lobes)

# If the PYTHON_LIBRARY_DIR is not specified, install in lib dir
if(NOT DEFINED ${PYTHON_LIBRARY_DIR}) 
   set(PYTHON_LIBRARY_DIR ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)
endif()

# Install pyeverybeam in site-packages directory
install(TARGETS pyeverybeam
  COMPONENT python
  LIBRARY DESTINATION ${PYTHON_LIBRARY_DIR}
)

# ${PYTHON_LIBRARY_DIR}