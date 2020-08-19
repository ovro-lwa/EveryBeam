#!/bin/sh

export PATH=$EXTRA_PATH:$PATH
python3 -B `dirname "${0}"`/compare_stationresponse.py

