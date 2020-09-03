#!/bin/sh

export PATH=$EXTRA_PATH:$PATH
python3 -B `dirname "${0}"`/run_oskar_simulation.py
add_beaminfo.py testdata.ms ${DATA_DIR}/skalowmini-coef.tm
