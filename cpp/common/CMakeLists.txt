# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

install(FILES constants.h casautils.h mathutils.h types.h fftresampler.h
              sphericalharmonics.h
        DESTINATION "include/${CMAKE_PROJECT_NAME}/common")
