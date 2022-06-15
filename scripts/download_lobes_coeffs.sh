#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading and extracting the LOBES coefficients into
# the current directory. When invoked from within CMake, CMake should set
# the proper WORKING_DIRECTORY.
set -e

# Download the h5 coefficient files in case they're not present (-nc option) in the coeffs/lobes directory
wget -q -r -nc -nH --no-parent --cut-dirs=4 --accept-regex=LOBES_.*\.h5 -e robots=off https://support.astron.nl/software/lobes/
