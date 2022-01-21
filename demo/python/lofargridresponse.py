# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import everybeam as eb
import numpy as np

"""
Demo script illustrating how the EveryBeam Python
bindings can be used to compute the beam response on
a grid for LOFAR observations.

Script requires an LOFAR Measurement Set in order to run.
A minimal LOFAR MSet can be downloaded here:
https://www.astron.nl/citt/EveryBeam/lba.MS.tar.bz2
"""

# Set path to LOFAR LBA MS and load telescope
ms_path = os.path.join(os.environ["DATA_DIR"], "LOFAR_LBA_MOCK.ms")

telescope = eb.load_telescope(ms_path)
assert type(telescope) == eb.LOFAR

# Set time and freq for which beam should be evaluated
time = 4.92183348e09
freq = 57884216.30859375
station_id = 23

# Specify the settings of the grid (/image)
gs = eb.GridSettings()
gs.width = gs.height = 4
gs.ra = -1.44194878
gs.dec = 0.85078091
gs.dl = gs.dm = 0.5 * np.pi / 180.0
gs.l_shift = gs.m_shift = 0.0

# Get the gridded response for all stations at once
grid_response_all = telescope.gridded_response(gs, time, freq)

# Check whether the returned numpy array has the expected shape
assert grid_response_all.shape == (telescope.nr_stations, gs.width, gs.height, 2, 2)

# Get the response for a specific station (station_id)
grid_response = telescope.gridded_response(gs, time, freq, station_id)

# Returned response should of course match corresponding entry
# in "grid_response_all"
np.testing.assert_allclose(grid_response, grid_response_all[station_id, ...], rtol=1e-6)

# Also, at the center pixel (2,2) the retrieved response should be reproduced
# by a call to station_response at the phase centre ra and dec
station_response = telescope.station_response(time, station_id, freq, gs.ra, gs.dec)
np.testing.assert_allclose(grid_response[2, 2, ...], station_response, rtol=1e-6)
