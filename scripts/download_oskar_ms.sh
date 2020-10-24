#!/bin/sh
#
# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading a mock OSKAR Measurement set

set -e

SCRIPT_PATH=$(dirname "$0")
cd $SCRIPT_PATH

# Move up to parent folder which contains the source
cd ..
mkdir -p test_data
cd test_data/

# TODO: complete name of tarfile to be downloaded
OSKAR_MOCK_ARCHIVE=OSKAR_MOCK.tar.bz2
OSKAR_MOCK_MS=OSKAR_MOCK.ms

if [ ! -f ${OSKAR_MOCK_MS}/table.f0 ] ; then

    if [ ! -f "$OSKAR_MOCK_ARCHIVE" ]; then
	wget -q https://www.astron.nl/citt/EveryBeam/OSKAR-single-timeslot.tar.bz2 -O $OSKAR_MOCK_ARCHIVE
    fi
    
    mkdir -p $OSKAR_MOCK_MS
    
    tar -xf $OSKAR_MOCK_ARCHIVE  -C $OSKAR_MOCK_MS --strip-components=1
    rm $OSKAR_MOCK_ARCHIVE
fi
