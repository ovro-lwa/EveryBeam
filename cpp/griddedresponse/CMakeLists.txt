# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

install(
  FILES griddedresponse.h
        lofargrid.h
        aartfaacgrid.h
        oskargrid.h
        dishgrid.h
        mwagrid.h
        phasedarraygrid.h
        skamidgrid.h
  DESTINATION "include/${CMAKE_PROJECT_NAME}/griddedresponse")
