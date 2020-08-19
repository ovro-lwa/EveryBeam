import os
import sys
import numpy as np
from read_oskar_beams import read_oskar_beams
import run_oskar_simulation
import subprocess


tolerance = float(os.environ["TOLERANCE"]) if "TOLERANCE" in os.environ else 0.0
npixels = int(os.environ["NPIXELS"]) if "NPIXELS" in os.environ else 256


run_oskar_simulation.main(npixels)
subprocess.check_call(["add_beaminfo.py", "skalowmini-coef.MS", "skalowmini-coef.tm"])
subprocess.check_call(["oskar_csv_to_hdf5.py", "skalowmini-coef.tm", "oskar.h5"])
subprocess.check_call(["./make_station_response_image", str(npixels)])

A = read_oskar_beams()
B = np.load('station-response.npy')

if tolerance:
    difference = np.nanmax(np.abs(A - B))
    if difference > tolerance:
        sys.exit(
            "Difference between OSKAR and EveryBeam spherical wave model is {}, which is larger than the tolerance {}".format(
                difference, tolerance
            )
        )

