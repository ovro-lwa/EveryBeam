# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import numpy as np
from read_oskar_beams import read_oskar_beams
from utils import check_tolerance
import run_oskar_simulation
import subprocess


tolerance = float(os.environ["TOLERANCE"]) if "TOLERANCE" in os.environ else 0.0
npixels = int(os.environ["NPIXELS"]) if "NPIXELS" in os.environ else 256
data_dir = os.environ["DATA_DIR"]


run_oskar_simulation.main(npixels)
subprocess.check_call(["add_beaminfo.py", "skalowmini-coef.MS", os.path.join(data_dir, "skalowmini-coef.tm")])
subprocess.check_call(["oskar_csv_to_hdf5.py", os.path.join(data_dir, "skalowmini-coef.tm"), "oskar.h5"])
subprocess.check_call(["./make_station_response_image", str(npixels)])

A = read_oskar_beams()
B = np.load('station-response.npy')

check_tolerance(tolerance, A, B)
