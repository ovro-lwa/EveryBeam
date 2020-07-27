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

VLA_MOCK_ARCHIVE=VLA_ARCHIVE.tar.bz2
VLA_MOCK_MS=VLA_MOCK.ms

if [ ! -f "$VLA_MOCK_ARCHIVE" ]; then
    wget https://www.astron.nl/citt/EveryBeam/small-vla-set.tar.bz2 -O $VLA_MOCK_ARCHIVE
fi

if [ -d $VLA_MOCK_MS ]
then
    echo "Directory already exists"
else
    mkdir $VLA_MOCK_MS
fi

tar -xf $VLA_MOCK_ARCHIVE  -C $VLA_MOCK_MS --strip-components=1
rm $VLA_MOCK_ARCHIVE