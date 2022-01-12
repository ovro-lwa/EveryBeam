# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# NOTE: the SKA-MID responses are compiled into a separate library, much like the oskar shared library
# for SKA-LOW. This simplifies a potential refactoring of the SKA telescope responses into a stand-alone library

add_library(skamidbeam SHARED skamidanalyticalresponse.cc)

string(TOLOWER ${CMAKE_PROJECT_NAME} projectname)
set_target_properties(skamidbeam PROPERTIES LIBRARY_OUTPUT_NAME
                                            "${projectname}-skamidbeam")
target_include_directories(skamidbeam PRIVATE ${Boost_INCLUDE_DIRS})

# When other targets within this project link against the skamid target, make sure
# they can find the include files.
target_include_directories(
  skamidbeam PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# install libeverybeam-skamidbeam.so
install(
  TARGETS skamidbeam
  EXPORT EveryBeamTargets
  DESTINATION lib)
