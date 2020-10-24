#!/bin/sh
#
# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading a mock LOFAR Measurement set

set -e

SCRIPT_PATH=$(dirname "$0")
cd $SCRIPT_PATH

# Move up to parent folder which contains the source
cd ..
mkdir -p test_data
cd test_data/

LOFAR_MOCK_ARCHIVE=LOFAR_HBA_ARCHIVE.tar.bz2
LOFAR_HBA_MOCK_MS=LOFAR_HBA_MOCK.ms

if [ ! -f ${LOFAR_HBA_MOCK_MS}/table.f0 ] ; then

    if [ ! -f "$LOFAR_MOCK_ARCHIVE" ]; then
	wget -q https://www.astron.nl/citt/EveryBeam/L258627-one-timestep.tar.bz2 -O $LOFAR_MOCK_ARCHIVE
    fi

    mkdir -p $LOFAR_HBA_MOCK_MS

    tar -xf $LOFAR_MOCK_ARCHIVE  -C $LOFAR_HBA_MOCK_MS --strip-components=1
    rm $LOFAR_MOCK_ARCHIVE
fi
