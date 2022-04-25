#!/bin/sh
# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script for downloading a measurement set used for KL screen fitting.


set -e

# Download SKA-MID mset
SCREEN_FITTING_MS_ARCHIVE=SCREEN_FITTING_MS.tar.bz2
SCREEN_FITTING_MS=SCREEN_FITTING.ms

if [ ! -f ${SCREEN_FITTING_MS}/table.f0 ] ; then

    if [ ! -f "$SCREEN_FITTING_MS_ARCHIVE" ]; then
	wget -q https://www.astron.nl/citt/ci_data/EveryBeam/screentest_ms.tar.gz -O $SCREEN_FITTING_MS_ARCHIVE
    fi

    rm -rf $SCREEN_FITTING_MS
    tar xf $SCREEN_FITTING_MS_ARCHIVE
    mv screen_test.ms $SCREEN_FITTING_MS

fi
