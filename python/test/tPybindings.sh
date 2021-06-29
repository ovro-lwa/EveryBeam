#!/bin/bash
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Following variables are assumed to be specified in CMakeLists file
# - $LIB_DIR path to pyeverybeam shared library
# - $DATA_DIR path to the (mock) test data
# - $SCRIPTS_DIR path to the (download) scripts

export PYTHONPATH=$LIB_DIR

pytest -s --exitfirst ${DIR}/test_pybindings.py