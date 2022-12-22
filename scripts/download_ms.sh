#!/bin/sh
# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script for downloading a Measurement Set

set -e

ARCHIVE=$1
MS=$2

if [ ! -f $MS/table.f0 ] ; then
    if [ ! -f $ARCHIVE ]; then
        wget -q https://support.astron.nl/software/ci_data/EveryBeam/$ARCHIVE
    fi

    mkdir -p $MS
    tar -xf $ARCHIVE -C $MS --strip-components=1
    rm -f $ARCHIVE
fi
