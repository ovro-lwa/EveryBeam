#!/bin/sh
#
# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading a mock VLA Measurement set

set -e

SCRIPT_PATH=$(dirname "$0")
cd $SCRIPT_PATH

# Move up to parent folder which contains the source
cd ..
mkdir -p test_data
cd test_data/

MWA_COEFF_ARCHIVE=MWA_COEFF.tar.bz2
MWA_COEFF_H5=MWA_COEFF.H5

if [ ! -f "$MWA_COEFF_ARCHIVE" ]; then
    wget -q www.astron.nl/citt/EveryBeam/mwa_full_embedded_element_pattern.tar.bz2 -O $MWA_COEFF_ARCHIVE
fi

tar -xf $MWA_COEFF_ARCHIVE
rm $MWA_COEFF_ARCHIVE