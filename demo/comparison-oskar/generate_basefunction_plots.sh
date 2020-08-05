#!/bin/sh

export PATH=$EXTRA_PATH:$PATH
export TOLERANCE="1e-12"
python3 -B `dirname "${0}"`/generate_basefunction_plots.py

