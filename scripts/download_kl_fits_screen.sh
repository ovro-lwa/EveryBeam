#!/bin/sh
# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script for downloading a reference fits file containing Karhunen-Loève screens.

set -e
KL_SCREEN_FILE="kl_0.fits"
# Download Karhunen Loève screens in fits format
if [ ! -f "$KL_SCREEN_FILE" ]; then
    wget https://www.astron.nl/citt/EveryBeam/${KL_SCREEN_FILE}
fi