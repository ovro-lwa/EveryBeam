#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

export PATH=$EXTRA_PATH:$PATH
python3 -B `dirname "${0}"`/generate_basefunction_plots.py
