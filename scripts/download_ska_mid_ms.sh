#!/bin/sh
# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script for downloading a minimal SKA-MID measurement set.
# The original simulated SKA-MID measurement set can be retrieved from Google Cloud Storage:
#
# gsutil -m rsync -r gs://ska1-simulation-data/simulations/continuum_simulations_SP-1331/mid/SKA_MID_SIM_custom_B2_dec_-45.0_nominal_nchan100.ms .
#
# This measurement set was reduced via two DP3 calls:
# DP3 msin=SKA_MID_SIM_custom_B2_dec_-45.0_nominal_nchan100.ms 'filter.startchan=0' filter.nchan=1 msout=SKA_MID_SIM_MOCK_TMP.ms steps=[filter] filter.remove=true msout.overwrite=true
# DP3 msin=SKA_MID_SIM_MOCK_TMP.ms msout=SKA_MID_SIM_MOCK.ms msin.ntimes=1 steps=[]


set -e

# Download SKA-MID mset
SKA_MID_MOCK_ARCHIVE=SKA_MID_MOCK.tar.bz2
SKA_MID_MOCK_MS=SKA_MID_MOCK.ms

if [ ! -f ${SKA_MID_MOCK_MS}/table.f0 ] ; then

    if [ ! -f "$SKA_MID_MOCK_ARCHIVE" ]; then
        wget -q https://support.astron.nl/software/ci_data/EveryBeam/SKA_MID_SIM-1channel-1timestep.tar.bz2 -O $SKA_MID_MOCK_ARCHIVE
    fi

    mkdir -p $SKA_MID_MOCK_MS

    tar -xf $SKA_MID_MOCK_ARCHIVE  -C $SKA_MID_MOCK_MS --strip-components=1
    rm -f $SKA_MID_MOCK_ARCHIVE
fi
