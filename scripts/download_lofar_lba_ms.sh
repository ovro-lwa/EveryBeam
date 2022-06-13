#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading a mock LOFAR (LBA) Measurement set

set -e

LOFAR_MOCK_ARCHIVE=LOFAR_LBA_ARCHIVE.tar.bz2
LOFAR_LBA_MOCK_MS=LOFAR_LBA_MOCK.ms

if [ ! -f ${LOFAR_LBA_MOCK_MS}/table.f0 ] ; then

    if [ ! -f "$LOFAR_MOCK_ARCHIVE" ]; then
        wget -q https://support.astron.nl/software/ci_data/EveryBeam/lba.MS.tar.bz2 -O $LOFAR_MOCK_ARCHIVE
    fi

    mkdir -p $LOFAR_LBA_MOCK_MS

    tar -xf $LOFAR_MOCK_ARCHIVE  -C $LOFAR_LBA_MOCK_MS --strip-components=1
    rm -f $LOFAR_MOCK_ARCHIVE
fi
