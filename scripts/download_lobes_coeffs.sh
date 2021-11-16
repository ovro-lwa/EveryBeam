#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading and extracting the LOBES coefficients
# When invoked from within CMAKE, the WORKING_DIRECTORY is
# assumed to be the CMAKE_BINARY_DIR
set -e

# Make directory coeffs/lobes in source directory
mkdir -p coeffs/lobes
cd coeffs/lobes

# Download the h5 coefficient files in case they're not present (-nc option) in the coeffs/lobes directory
wget -q -r -nc -nH --no-check-certificate --no-parent --cut-dirs=3 --accept-regex=LOBES_.*\.h5 --reject index.html -e robots=off https://www.astron.nl/citt/EveryBeam/lobes/
