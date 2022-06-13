#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading a mock VLA Measurement set

set -e

MWA_MOCK_ARCHIVE=MWA_ARCHIVE.tar.bz2
MWA_MOCK_MS=MWA_MOCK.ms

if [ ! -f ${MWA_MOCK_MS}/table.f1 ]; then

    if [ ! -f "$MWA_MOCK_ARCHIVE" ]; then
        wget -q https://support.astron.nl/software/ci_data/EveryBeam/MWA-single-timeslot.tar.bz2 -O $MWA_MOCK_ARCHIVE
    fi

    mkdir -p $MWA_MOCK_MS

    tar -xf $MWA_MOCK_ARCHIVE  -C $MWA_MOCK_MS --strip-components=1
    rm -f $MWA_MOCK_ARCHIVE
fi
