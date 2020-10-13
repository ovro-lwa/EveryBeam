#!/bin/sh
#
# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading a mock LOFAR (LBA) Measurement set

set -e

SCRIPT_PATH=$(dirname "$0")
cd $SCRIPT_PATH

# Move up to parent folder which contains the source
cd ..
mkdir -p test_data
cd test_data/

LOFAR_MOCK_ARCHIVE=LOFAR_LBA_ARCHIVE.tar.bz2
LOFAR_LBA_MOCK_MS=LOFAR_LBA_MOCK.ms

if [ ! -f "$LOFAR_MOCK_ARCHIVE" ]; then
    wget -q https://www.astron.nl/citt/EveryBeam/lba.MS.tar.bz2 -O $LOFAR_MOCK_ARCHIVE
fi

if [ -d $LOFAR_LBA_MOCK_MS ]
then
    echo "Directory already exists"
else
    mkdir $LOFAR_LBA_MOCK_MS
fi

tar -xf $LOFAR_MOCK_ARCHIVE  -C $LOFAR_LBA_MOCK_MS --strip-components=1
rm $LOFAR_MOCK_ARCHIVE