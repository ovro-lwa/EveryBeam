#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

set -e

SRC_DIR="$( cd "$( dirname "${0}" )" >/dev/null 2>&1 && pwd )"

atermoptions="-aterm-config ${SRC_DIR}/aterm.conf" #  -save-aterms"

intervaloptions="-interval 0 100"

# Workaround to generate beam image oskar-1-1-beam-I.fits
wsclean ${intervaloptions} -name oskar-1 -no-dirty -scale 2amin -size 300 300 ${atermoptions} -use-idg -niter 1 ${MS}

# run wsclean on OSKAR generated data, without and with beam correction
wsclean ${intervaloptions} -name oskar-0 -data-column DATA -pol IQUV -no-dirty -scale 2amin -size 300 300 ${MS}
wsclean ${intervaloptions} -name oskar-1 -data-column DATA -pol IQUV -no-dirty -scale 2amin -size 300 300 ${atermoptions} -use-idg ${MS}

# use wsclean to generate model data without beam, and to image it without beam correction
wsclean ${intervaloptions} -predict -pol IQUV -use-idg ${MS}
wsclean ${intervaloptions} -name oskar-0-0 -data-column MODEL_DATA -pol IQUV -no-dirty -scale 2amin -size 300 300 -use-idg ${MS}
# wsclean ${intervaloptions} -name beam-0-1 -data-column MODEL_DATA -pol IQUV -no-dirty -scale 2amin -size 300 300 ${atermoptions} -use-idg ${MS}

# Workaround to generate beam image oskar-1-1-beam-I.fits
wsclean ${intervaloptions} -name oskar-1-1 -no-dirty -scale 2amin -size 300 300 ${atermoptions} -use-idg -niter 1 ${MS}

# use wsclean to generate model data with beam, and to image it with and without beam correction
wsclean ${intervaloptions} -predict -pol IQUV ${atermoptions} -use-idg ${MS}
wsclean ${intervaloptions} -name oskar-1-1 -data-column MODEL_DATA -pol IQUV -no-dirty -scale 2amin -size 300 300 ${atermoptions} -use-idg ${MS}
wsclean ${intervaloptions} -name oskar-1-0 -data-column MODEL_DATA -pol IQUV -no-dirty -scale 2amin -size 300 300 -use-idg ${MS}
