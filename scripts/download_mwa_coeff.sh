#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading the MWA coefficients

set -e

SCRIPT_PATH=$(dirname "$0")
cd $SCRIPT_PATH

# Move up to parent folder which contains the source
cd ..
mkdir -p test_data
cd test_data/

MWA_COEFF_ARCHIVE=MWA_COEFF.tar.bz2
MWA_COEFF_H5=mwa_full_embedded_element_pattern.h5

if [ ! -f ${MWA_COEFF_H5} ] ; then

    wget -q www.astron.nl/citt/EveryBeam/mwa_full_embedded_element_pattern.tar.bz2 -O ${MWA_COEFF_ARCHIVE}

    tar -xjf $MWA_COEFF_ARCHIVE
    rm $MWA_COEFF_ARCHIVE
fi
