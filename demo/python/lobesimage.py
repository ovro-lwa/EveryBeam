# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

from everybeam import load_telescope, LOFAR, Options, thetaphi2cart
import matplotlib.pyplot as plt
import numpy as np
import logging

"""
Demo script illustrating how the EveryBeam Python
bindings can be used to generate all-sky images
for the beam response for a LOFAR observation.

Script requires a LOFAR Measurement Set in order to run.
A minimal LOFAR MSet can be downloaded here:
https://www.astron.nl/citt/EveryBeam/lba.MS.tar.bz2
"""

ms_path = os.path.join(os.environ["DATA_DIR"], "LOFAR_LBA_MOCK.ms")

# Response settings
mode = "station"  # Change to "element" to obtain the element response
time = 4.92183348e09
frequency = 57812500.0
station_id = 20
element_id = 0

# Same station0 direction as in python/test
station0 = np.array([0.655743, -0.0670973, 0.751996])
is_local = True  # use local coords.
rotate = True
response_model = "lobes"

# Load telescope
telescope = load_telescope(ms_path, element_response_model=response_model)
station_name = telescope.station_name(station_id)

if station_name not in ["CS302LBA", "SE607LBA"] and response_model == "lobes":
    logging.warning(
        f"LOBEs response requested, but not available for station {station_name}. EveryBeam defaults back to Hamaker"
    )

# Make coordinates for mesh
x_v = np.linspace(-1.0, 1.0, 128)
y_v = np.linspace(-1.0, 1.0, 128)

response = np.empty((y_v.size, x_v.size, 2, 2), dtype=np.cdouble)
response.fill(np.nan)
for i, x in enumerate(x_v):
    for j, y in enumerate(y_v):
        if (x ** 2 + y ** 2) <= 1.0:
            # Compute theta/phi and resulting direction vector
            theta = np.arcsin(np.sqrt(x * x + y * y))
            phi = np.arctan2(y, x)
            direction = thetaphi2cart(theta, phi)
            if mode == "element":
                response[j, i, :, :] = telescope.element_response(
                    time,
                    station_id,
                    element_id,
                    frequency,
                    direction,
                    is_local,
                    rotate=rotate,
                )
            elif mode == "station":
                response[j, i, :, :] = telescope.station_response(
                    time, station_id, frequency, direction, station0, rotate=rotate
                )
            else:
                raise Exception(
                    "Unrecognized response mode. Must be either station or element"
                )

# Plotting
X, Y = np.meshgrid(x_v, y_v)
label = ("X", "Y")
fig, axs = plt.subplots(2, 4, figsize=(10, 5), sharex=True, sharey=True)
fig.suptitle(f"{mode} response for station {station_name}: {response_model}")
for i, row in enumerate(range(2)):
    for j, col in enumerate(range(2)):
        im1 = axs[row, col].pcolor(X, Y, np.abs(response[:, :, row, col]))
        axs[row, col].set_aspect("equal", "box")
        axs[row, col].set_title(f"abs({label[j]}{label[i]})")
        fig.colorbar(im1, ax=axs[row, col])

        im2 = axs[row, col + 2].pcolor(X, Y, np.angle(response[:, :, row, col]))
        axs[row, col + 2].set_aspect("equal", "box")
        axs[row, col + 2].set_title(f"angle({label[j]}{label[i]})")
        fig.colorbar(im2, ax=axs[row, col + 2])

plt.show()
