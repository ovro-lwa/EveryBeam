# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import numpy as np
import matplotlib.pyplot as plt
import run_oskar
from oskar_utils import check_tolerance, generate_oskar_csv, read_oskar_beams
import subprocess

# Check and set parameters
l_max = int(os.environ["MAX_ORDER"]) if "MAX_ORDER" in os.environ else 2
apply_transpose = (
    os.environ["APPLY_TRANSPOSE"].upper() in ("1", "ON", "TRUE")
    if "APPLY_TRANSPOSE" in os.environ
    else False
)
tolerance = float(os.environ["TOLERANCE"]) if "TOLERANCE" in os.environ else 0.0
npixels = int(os.environ["NPIXELS"]) if "NPIXELS" in os.environ else 256

plt.figure(figsize=(10, 6))

# Freqs (in MHz)
freqs = [50]
output_dir = "telescope.tm"
element_id = 0

root_name = f"basefunctions_{freqs[0]:03}"
for em_idx in range(2):
    for basefunction_idx in range(l_max * (l_max + 2)):
        plt.clf()

        generate_oskar_csv(basefunction_idx, freqs, output_dir, element_id, em_idx)
        run_oskar.main(npixels)
        A = read_oskar_beams(root_name)

        plt.subplot(2, 4, 1)
        plt.imshow(np.abs(A[:, :, 0, 0]).T, clim=(0, 0.25), origin="lower")
        plt.colorbar()
        plt.title("abs(Etheta)")
        plt.ylabel("OSKAR")

        plt.subplot(2, 4, 2)
        plt.imshow(np.abs(A[:, :, 0, 1]).T, clim=(0, 0.25), origin="lower")
        plt.colorbar()
        plt.title("abs(Ephi)")

        plt.subplot(2, 4, 3)
        plt.imshow(
            np.angle(A[:, :, 0, 0]).T,
            clim=(-np.pi, np.pi),
            cmap="twilight",
            origin="lower",
        )
        plt.colorbar()
        plt.title("angle(Etheta)")

        plt.subplot(2, 4, 4)
        plt.imshow(
            np.angle(A[:, :, 0, 1]).T,
            clim=(-np.pi, np.pi),
            cmap="twilight",
            origin="lower",
        )
        plt.colorbar()
        plt.title("angle(Ephi)")

        l = int(np.sqrt(basefunction_idx + 1))
        m = basefunction_idx - l * l + 1 - l
        s = em_idx

        if not apply_transpose:
            generate_oskar_csv(basefunction_idx, freqs, output_dir, element_id, em_idx)
        else:
            # flip the sign of m
            generate_oskar_csv(l * l - 1 + l - m, freqs, output_dir, element_id, em_idx)

        subprocess.check_call(["oskar_csv_to_hdf5.py", "telescope.tm", "oskar.h5"])
        subprocess.check_call(["make_element_response_image", str(npixels)])

        B = np.load("response.npy")

        if apply_transpose:
            B *= (-1) ** (m + 1)

        plt.subplot(2, 4, 5)
        plt.imshow(np.abs(B[:, :, 0, 0]).T, clim=(0, 0.25), origin="lower")
        plt.colorbar()
        plt.title("abs(Etheta)")
        plt.ylabel("EveryBeam")

        plt.subplot(2, 4, 6)
        plt.imshow(np.abs(B[:, :, 0, 1]).T, clim=(0, 0.25), origin="lower")
        plt.colorbar()
        plt.title("abs(Ephi)")

        plt.subplot(2, 4, 7)
        plt.imshow(
            np.angle(B[:, :, 0, 0]).T,
            clim=(-np.pi, np.pi),
            cmap="twilight",
            origin="lower",
        )
        plt.colorbar()
        plt.title("angle(Etheta)")

        plt.subplot(2, 4, 8)
        plt.imshow(
            np.angle(B[:, :, 0, 1]).T,
            clim=(-np.pi, np.pi),
            cmap="twilight",
            origin="lower",
        )
        plt.colorbar()
        plt.title("angle(Ephi)")
        plt.gcf().suptitle("l = {}, m = {}, s = {}".format(l, m, s))
        plt.savefig("basefunction{}-{}".format(basefunction_idx, em_idx))

        check_tolerance(tolerance, A, B)
