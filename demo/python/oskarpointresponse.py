# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import everybeam as eb
import numpy as np

"""
Demo script illustrating how the EveryBeam Python
bindings can be used to compute beam responses.

Script requires an OSKAR Measurement Set in order to run.
A minimal OSKAR MSet can be downloaded here:
https://www.astron.nl/citt/EveryBeam/OSKAR-single-timeslot.tar.bz2
"""

# Set path to (OSKAR) MS
ms_path = "../../test_data/OSKAR_MOCK.ms"

# Use differential beam?
use_differential_beam = False

# Set element response model
element_response_model = "skala40_wave"

# Load the telescope
telescope = eb.load_telescope(
    ms_path,
    use_differential_beam=use_differential_beam,
    element_response_model=element_response_model,
)

# Is this an OSKAR telescope?
assert type(telescope) == eb.OSKAR

# Set properties for the beam evaluation
time = 4.45353e09
station_id = 0
freq = 5.0e07

# Point of interest (given in ITRF coordinates)
dir_itrf = np.array([-0.203137, 0.841818, -0.500078])

# Array factor for station 0, since we are evaluating at the phase centre (delay direction and
# direction of interest coincide), a unity matrix should be returned
array_factor_phase_centre = telescope.array_factor(
    time, station_id, freq, dir_itrf, dir_itrf
)
np.testing.assert_allclose(
    array_factor_phase_centre, np.eye(2, dtype=np.complex64), rtol=1e-8
)

# Just a rather arbitrarily chosen delay direction
delay_dir_itrf = dir_itrf + np.array([-0.1, 0.5, 0.025])
array_factor = telescope.array_factor(time, station_id, freq, dir_itrf, delay_dir_itrf)

# Element beam for station 0
element_response = telescope.element_response(time, station_id, freq, dir_itrf)

# Full beam for station 0
response = telescope.station_response(time, station_id, freq, dir_itrf, delay_dir_itrf)

# Full beam should match product of array_factor and element_response
np.testing.assert_allclose(
    np.matmul(array_factor, element_response), response, rtol=1e-6
)
