#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

export PATH=$EXTRA_PATH:$PATH
python3 -B `dirname "${0}"`/run_oskar_simulation.py
add_beaminfo.py testdata.ms ${DATA_DIR}/skalowmini-coef.tm
