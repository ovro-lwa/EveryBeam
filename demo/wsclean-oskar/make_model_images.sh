#!/bin/sh

export PATH=$EXTRA_PATH:$PATH
python3 -B `dirname "${0}"`/make_model_images.py

