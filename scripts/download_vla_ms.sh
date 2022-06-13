#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading a mock VLA Measurement set

set -e

VLA_MOCK_ARCHIVE=VLA_ARCHIVE.tar.bz2
VLA_MOCK_MS=VLA_MOCK.ms

if [ ! -f ${VLA_MOCK_MS}/table.f1 ] ; then

    if [ ! -f "$VLA_MOCK_ARCHIVE" ]; then
        wget -q https://support.astron.nl/software/ci_data/EveryBeam/small-vla-set.tar.bz2 -O $VLA_MOCK_ARCHIVE
    fi

    mkdir -p $VLA_MOCK_MS

    tar -xf $VLA_MOCK_ARCHIVE  -C $VLA_MOCK_MS --strip-components=1
    rm -f $VLA_MOCK_ARCHIVE
fi
