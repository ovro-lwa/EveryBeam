# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

project(lobes_pybind11)

# Create the binding library
pybind11_add_module(pylobes 
                    lobes.cc
                    )


target_include_directories(pylobes PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(pylobes PUBLIC ${HDF5_LIBRARIES})
