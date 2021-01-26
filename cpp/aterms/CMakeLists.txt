# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

install (FILES
atermconfig.h
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
DESTINATION "include/${CMAKE_PROJECT_NAME}/aterms")