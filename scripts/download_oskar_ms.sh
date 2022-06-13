#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading a mock OSKAR Measurement set

set -e

# TODO: complete name of tarfile to be downloaded
OSKAR_MOCK_ARCHIVE=OSKAR_MOCK.tar.bz2
OSKAR_MOCK_MS=OSKAR_MOCK.ms

if [ ! -f ${OSKAR_MOCK_MS}/table.f0 ] ; then

    if [ ! -f "$OSKAR_MOCK_ARCHIVE" ]; then
        wget -q https://support.astron.nl/software/ci_data/EveryBeam/OSKAR-single-timeslot.tar.bz2 -O $OSKAR_MOCK_ARCHIVE
    fi

    mkdir -p $OSKAR_MOCK_MS

    tar -xf $OSKAR_MOCK_ARCHIVE  -C $OSKAR_MOCK_MS --strip-components=1
    rm -f $OSKAR_MOCK_ARCHIVE
fi
