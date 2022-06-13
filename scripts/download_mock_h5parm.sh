#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

# Script for downloading a mock h5parm file
# for testing the H5ParmATerms

set -e

MOCK_H5PARM=MOCK_H5PARM.h5
if [ ! -f ${MOCK_H5PARM} ] ; then
    wget -q https://support.astron.nl/software/ci_data/EveryBeam/${MOCK_H5PARM}
fi