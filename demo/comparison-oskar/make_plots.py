# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
from matplotlib import pyplot as plt
from read_oskar_beams import read_oskar_beams

A = read_oskar_beams()
B = np.load('station-response.npy')

B = B-A

#B = np.load('response.npy')

for ant_xy_idx in range(2):

    plt.subplot(2, 4, 1 + ant_xy_idx*4)
    plt.imshow(np.abs(B[:, :, ant_xy_idx, 0]).T, origin="lower")
    plt.colorbar()
    plt.title("abs(Etheta)")

    plt.subplot(2, 4, 2 + ant_xy_idx*4)
    plt.imshow(np.abs(B[:, :, ant_xy_idx, 1]).T, origin="lower")
    plt.colorbar()
    plt.title("abs(Ephi)")

    plt.subplot(2, 4, 3 + ant_xy_idx*4)
    plt.imshow(
        np.angle(B[:, :, ant_xy_idx, 0]).T, clim=(-np.pi, np.pi), cmap="twilight", origin="lower"
    )
    plt.colorbar()
    plt.title("angle(Etheta)")

    plt.subplot(2, 4, 4 + ant_xy_idx*4)
    plt.imshow(
        np.angle(B[:, :, ant_xy_idx, 1]).T, clim=(-np.pi, np.pi), cmap="twilight", origin="lower"
    )
    plt.colorbar()
    plt.title("angle(Ephi)")

plt.show()
