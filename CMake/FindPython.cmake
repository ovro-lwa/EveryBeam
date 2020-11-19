# - Try to find the Python interpreter, Python header files and libraries.
# This macro effectively wraps the FindPythonInterp and FindPythonLibs macros
# provided by CMake.
#
# In addition to the variables that are set by FindPythonInterp and
# FindPythonLibs, this will define:
#  PYTHON_FOUND        - system has Python interpreter, Python headers
#                        files and libraries
#  PYTHON_INCLUDE_DIRS - path to the Python header files
#  PYTHON_BUILD_DIR    - build directory for Python extensions (cached)
#  PYTHON_INSTALL_DIR  - installation directory for Python extensions (cached)

# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Set options string to pass to the find_package() commands below.
set(_options ${Python_FIND_VERSION})
if(Python_FIND_VERSION_EXACT)
  list(APPEND _options EXACT)
endif(Python_FIND_VERSION_EXACT)
if(Python_FIND_QUIETLY)
  list(APPEND _options QUIET)
endif(Python_FIND_QUIETLY)
if(Python_FIND_REQUIRED)
  list(APPEND _options REQUIRED)
endif(Python_FIND_REQUIRED)

# Search for the Python interpreter.
find_package(PythonInterp ${_options})

# Search for the Python header files and libraries.
find_package(PythonLibs ${_options})

# Set PYTHON_INCLUDE_DIRS variable, because FindPythonLibs does not do it.
set(PYTHON_INCLUDE_DIRS "${PYTHON_INCLUDE_PATH}")

# PythonInstall sets PYTHON_BUILD_DIR and PYTHON_INSTALL_DIR
include(PythonInstall)

# Set PYTHON_FOUND to TRUE if both Python interpreter and libraries are found.
set(PYTHON_FOUND FALSE)
if(PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND)
  set(PYTHON_FOUND TRUE)
endif(PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND)
