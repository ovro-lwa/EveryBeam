# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

install(
  FILES pointresponse.h
        phasedarraypoint.h
        aartfaacpoint.h
        lofarpoint.h
        oskarpoint.h
        dishpoint.h
        mwapoint.h
  DESTINATION "include/${CMAKE_PROJECT_NAME}/pointresponse")
