#!/bin/sh

export PATH=$EXTRA_PATH:$PATH
python3 -B `dirname "${0}"`/generate_basefunction_plots.py

