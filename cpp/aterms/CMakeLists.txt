# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(
  aterms OBJECT
  atermconfig.cc
  atermbase.cc
  atermresampler.cc
  everybeamaterm.cc
  fitsatermbase.cc
  fitsaterm.cc
  dldmaterm.cc
  pafbeamterm.cc
  h5parmaterm.cc
  klfitter.cc
  klfittingaterm.cc
  fourierfitter.cc
  fourierfittingaterm.cc
  cache.cc)

# The aterms object library is included in the SHARED library everybeam
# Compilation must therefore generate relocatable code
set_property(TARGET aterms PROPERTY POSITION_INDEPENDENT_CODE ON)

set_property(TARGET aterms PROPERTY CXX_STANDARD 17)

target_include_directories(
  aterms SYSTEM
  PRIVATE ${CMAKE_SOURCE_DIR}/external/xtensor/include
          ${CMAKE_SOURCE_DIR}/external/xtensor-blas/include
          ${CMAKE_SOURCE_DIR}/external/xtensor-fftw/include
          ${CMAKE_SOURCE_DIR}/external/xtl/include)

install(
  FILES atermconfig.h
        parsetprovider.h
        # May not need to expose these
        cache.h
        atermbase.h
        atermbeam.h
        atermresampler.h
        #
        fitsatermbase.h
        fitsaterm.h
        everybeamaterm.h
        dldmaterm.h
        pafbeamterm.h
        h5parmaterm.h
        klfittingaterm.h
        fourierfittingaterm.h
  DESTINATION "include/${CMAKE_PROJECT_NAME}/aterms")
