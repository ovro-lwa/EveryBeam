#!/bin/sh
# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script for downloading a reference h5 solutions file.
# This file is generated from the "screen_test.ms" dataset, using the following DP3 command:

# DP3 msin.datacolumn=DATA msout= msin=screen_test.ms steps=[ddecal] ddecal.type=ddecal ddecal.mode=scalarphase ddecal.h5parm=test_solutions.h5parm ddecal.maxiter=150 ddecal.nchan=10 ddecal.solint=1 ddecal.sourcedb=skymodel.txt ddecal.usebeammodel=True ddecal.beammode=array_factor ddecal.solveralgorithm=hybrid ddecal.stepsize=0.02
# Where the skymodel is the following:
# FORMAT = Name, Type, Patch, Ra, Dec, I, SpectralIndex='[]', LogarithmicSI, ReferenceFrequency='150616963.704427', MajorAxis, MinorAxis, Orientation
# # LSMTool history:
# # 2022-01-12 11:36:54: LOAD (from file '/beegfs/rafferty/Data/LOFAR/Screens/input_data_sim/sim_skymodel.txt')
# # 2022-01-12 11:36:54: SETPATCHPOSITIONS (method = 'mid')
#  , , Patch_0, 8:37:42.9518, 65.13.47.4993
#  , , Patch_1, 8:13:12.5713, 65.12.15.3123
#  , , Patch_10, 8:19:04.7255, 64.30.13.5853
#  , , Patch_2, 8:11:33.4617, 63.02.23.3212
#  , , Patch_3, 8:35:12.5497, 63.03.22.1235
#  , , Patch_4, 8:29:19.4111, 63.32.51.7385
#  , , Patch_8, 8:23:07.1265, 64.05.34.771
# s1, POINT, Patch_0, 8:37:42.9518, 65.13.47.4993, 5.0, [-0.076111092495773, -0.0419417198821972], false, 150616963.704427, 0, 0, 0
# s7, POINT, Patch_1, 8:13:12.5713, 65.12.15.3123, 5.0, [-0.133957875863411, -0.91057631218414], false, 150616963.704427, 0, 0, 0
# s4, POINT, Patch_10, 8:19:04.7255, 64.30.13.5853, 5.0, [-0.0368563521820315, -0.233815778288296], false, 150616963.704427, 0, 0, 0
# s5, POINT, Patch_2, 8:11:33.4617, 63.02.23.3212, 5.0, [-0.123957875863411, -0.89057631218414], false, 150616963.704427, 0, 0, 0
# s6, POINT, Patch_3, 8:35:12.5497, 63.03.22.1235, 5.0, [-0.066111092495773, -0.0429417198821972], false, 150616963.704427, 0, 0, 0
# s2, POINT, Patch_4, 8:29:19.4111, 63.32.51.7385, 5.0, [-0.0468212171772782, 0.116999569372316], false, 150616963.704427, 0, 0, 0
# s3, POINT, Patch_8, 8:23:07.1265, 64.05.34.771, 5.0, [-0.0303229488857636, -0.106658633984568], false, 150616963.704427, 0, 0, 0
# s3a, POINT, Patch_8, 8:22:07.1265, 64.05.34.771, 0.5, [-0.0303229488857636, -0.106658633984568], false, 150616963.704427, 0, 0, 0
# s3b, POINT, Patch_8, 8:24:07.1265, 64.05.34.771, 0.5, [-0.0468212171772782, 0.116999569372316], false, 150616963.704427, 0, 0, 0
# s3c, POINT, Patch_8, 8:23:07.1265, 64.25.34.771, 0.5, [-0.123957875863411, -0.89057631218414], false, 150616963.704427, 0, 0, 0
# s3d, POINT, Patch_8, 8:23:07.1265, 63.45.34.771, 0.5, [-0.0303229488857636, -0.106658633984568], false, 150616963.704427, 0, 0, 0


set -e
H5_SOLUTIONS_FILE="solutions.h5"
# Download Karhunen Loève screens in fits format
if [ ! -f "$H5_SOLUTIONS_FILE" ]; then
    wget https://support.astron.nl/software/ci_data/EveryBeam/${H5_SOLUTIONS_FILE}
fi