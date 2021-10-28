#!/usr/bin/env python3
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Generates all-sky zenith-centred beam patterns for SKALA-4 and EDA-2 antennas.
"""

import copy
import os
import subprocess

# pylint: disable=wrong-import-position
from oskar_utils import create_settings, get_telescope_settings


def main(npixels):
    """Main function."""
    # Name of the application to run, and a settings file for it.
    app1 = "oskar_sim_interferometer"
    app2 = "oskar_sim_beam_pattern"
    settings_path = "_temp_settings.ini"

    data_dir = os.environ["DATA_DIR"]

    # Define some basic observation parameters.
    ra0_deg = 0.0
    dec0_deg = -27.0
    length_sec = 60  # 3600.0*6
    # Define frequencies of interest (in MHz).
    freqs = [50]

    telescope_settings = get_telescope_settings(
        data_dir, ra0_deg, dec0_deg, length_sec, num_time_steps=60, swap_xy=False
    )

    # Copy the base settings dictionary.
    current_settings = copy.deepcopy(telescope_settings)

    # Loop over frequencies.
    for freq in freqs:
        # Create the settings file.
        settings = create_settings(app1, settings_path, current_settings)

        # Update output root path and frequency.
        settings["observation/start_frequency_hz"] = 1e6 * freq
        settings["interferometer/ms_filename"] = "skalowmini-coef.MS"

        # Run the app with the settings file.
        subprocess.call([app1, settings_path])

        # Create the settings file.
        settings = create_settings(app2, settings_path, current_settings)

        # Update frequency.
        root_path = "basefunctions" + f"_{int(freq):03}_MHz"
        settings["beam_pattern/root_path"] = root_path
        settings["observation/start_frequency_hz"] = 1e6 * freq
        settings["beam_pattern/coordinate_frame"] = "Horizon"
        settings["beam_pattern/beam_image/size"] = npixels
        settings["beam_pattern/station_outputs/fits_image/amp"] = True
        settings["beam_pattern/station_outputs/fits_image/phase"] = True
        settings["beam_pattern/station_outputs/fits_image/auto_power_real"] = False

        # Run the app with the settings file.
        subprocess.call([app2, settings_path])


if __name__ == "__main__":
    main(256)
