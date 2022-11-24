#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-only

# Script to install FFTW from source

set -euo pipefail

pushd /tmp

echo "Downloading & unpacking FFTW ${FFTW_VERSION}"
curl https://fftw.org/pub/fftw/fftw-${FFTW_VERSION}.tar.gz --output fftw-${FFTW_VERSION}.tar.gz

tar xf fftw-${FFTW_VERSION}.tar.gz
pushd fftw-${FFTW_VERSION}

echo "Configuring, building & installing FFTW ${FFTW_VERSION} to ${FFTW_DIR}"
./configure --prefix=/opt/fftw3 --enable-threads --enable-shared --enable-float
make -j $(nproc)
make install
popd

# Clean up to limit the size of the Docker image
echo "Cleaning up unnecessary files"
rm -r fftw-${FFTW_VERSION}
rm fftw-${FFTW_VERSION}.tar.gz
