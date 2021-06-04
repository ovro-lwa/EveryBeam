# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# TODO: probably, we only need/want to expose the telescope.h header
install(
  FILES dish.h
        lofar.h
        mwa.h
        lofar.h
        oskar.h
        phasedarray.h
        telescope.h
  DESTINATION "include/${CMAKE_PROJECT_NAME}/telescope")
