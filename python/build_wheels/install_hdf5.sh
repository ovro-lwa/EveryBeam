#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-only

# Script to install HDF5 from source

set -euo pipefail

echo "Installing zlib with yum"
yum -y install zlib-devel

pushd /tmp

# This seems to be needed to find libsz.so.2
ldconfig

echo "Downloading & unpacking HDF5 ${HDF5_VERSION}"
#                                   Remove trailing .*, to get e.g. '1.12' ↓
curl -fsSLO "https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_VERSION%.*}/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz"
tar -xzvf hdf5-${HDF5_VERSION}.tar.gz
pushd hdf5-${HDF5_VERSION}

echo "Configuring, building & installing HDF5 ${HDF5_VERSION} to ${HDF5_DIR}"
./configure --prefix ${HDF5_DIR} --enable-build-mode=production --with-szlib --enable-cxx
make -j $(nproc)
make install
popd

# Clean up to limit the size of the Docker image
echo "Cleaning up unnecessary files"
rm -r hdf5-${HDF5_VERSION}
rm hdf5-${HDF5_VERSION}.tar.gz

yum -y erase zlib-devel
