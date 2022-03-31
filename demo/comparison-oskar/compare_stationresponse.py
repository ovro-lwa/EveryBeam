# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Work-around issue https://github.com/astropy/astropy/issues/13007
from astropy.utils import iers
iers.conf.auto_download = False

import os
import numpy as np
from oskar_utils import check_tolerance, read_oskar_beams
import run_oskar_simulation
import subprocess
import everybeam as eb

"""
Compare the station response obtained with OSKAR to the station
response obtained via EveryBeam.
"""

tolerance = float(os.environ["TOLERANCE"]) if "TOLERANCE" in os.environ else 0.0
npixels = int(os.environ["NPIXELS"]) if "NPIXELS" in os.environ else 256
data_dir = os.environ["DATA_DIR"]
freq = 50e6
station_index = 0

# Obtain result with OSKAR
run_oskar_simulation.main(npixels)
subprocess.check_call(
    [
        "add_beaminfo.py",
        "skalowmini-coef.MS",
        os.path.join(data_dir, "skalowmini-coef.tm"),
    ]
)
subprocess.check_call(
    ["oskar_csv_to_hdf5.py", os.path.join(data_dir, "skalowmini-coef.tm"), "oskar.h5"]
)
# subprocess.check_call(["make_station_response_image", str(npixels)])

root_name = "basefunctions_050"
A = read_oskar_beams(root_name)


# Obtain results with EveryBeam

# Time and itrf vector of reference direction
# TODO: infer this from EveryBeam?
time = 4453522026
station0 = np.array([-0.4028537349, 0.7946940582, -0.4540597121])

# Antenna coordinate_system,
# TODO: infer these via EveryBeam?
# C++ counterpart is: station->GetAntenna()->coordinate_system_.axes.[p,q,r]
p = np.array([-0.8928654139, -0.4503236089, 0])
q = np.array([-0.2032142175, 0.4029167976, 0.89239119])
r = np.array([-0.4018648212, 0.7967852292, -0.451262633])

telescope = eb.load_telescope(
    "skalowmini-coef.MS", element_response_model="skala40_wave"
)
x_v = np.linspace(-1.0, 1.0, npixels)
y_v = np.linspace(-1.0, 1.0, npixels)
B = np.empty((y_v.size, x_v.size, 2, 2), dtype=np.cdouble)
B.fill(complex(np.nan, np.nan))

for i, x in enumerate(x_v):
    for j, y in enumerate(y_v):
        if (x ** 2 + y ** 2) <= 1.0:
            z = np.sqrt(1.0 - x ** 2 - y ** 2)
            direction = x * p + y * q + z * r
            B[i, j, :, :] = telescope.station_response(
                time, station_index, freq, direction, station0, rotate=False
            )

# Compare OSKAR (A) with EveryBeam (B)
check_tolerance(tolerance, A, B)
