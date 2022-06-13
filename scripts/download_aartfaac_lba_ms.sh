#!/bin/sh
# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Script for downloading a mock AARTFAAC measurement set

set -e

AARTFAAC_MOCK_ARCHIVE=LOFAR_LBA_ARCHIVE.tar.bz2
AARTFAAC_LBA_MOCK_MS=AARTFAAC_LBA_MOCK.ms

if [ ! -f ${AARTFAAC_LBA_MOCK_MS}/table.f0 ] ; then

    if [ ! -f "$LOFAR_MOCK_ARCHIVE" ]; then
        wget -q https://support.astron.nl/software/ci_data/EveryBeam/aartfaac.MS.tgz -O $AARTFAAC_MOCK_ARCHIVE
    fi

    mkdir -p $AARTFAAC_LBA_MOCK_MS

    tar -xf $AARTFAAC_MOCK_ARCHIVE  -C $AARTFAAC_LBA_MOCK_MS --strip-components=1
    rm -f $AARTFAAC_MOCK_ARCHIVE
fi
