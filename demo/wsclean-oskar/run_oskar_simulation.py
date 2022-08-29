#!/usr/bin/env python3
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Generates all-sky zenith-centred beam patterns for SKALA-4 and EDA-2 antennas.
"""

import copy
import os.path
import subprocess

import matplotlib

matplotlib.use("Agg")
# pylint: disable=wrong-import-position
import oskar
from oskar_utils import get_telescope_settings


def main():
    """
    Main function for running the OSKAR simulation.
    """

    # Name of the application to run, and a settings file for it.
    app1 = "oskar_sim_interferometer"
    settings_path = "_temp_settings.ini"
    data_dir = os.environ["DATA_DIR"]
    src_dir = os.path.dirname(__file__)

    # Define some basic observation parameters.
    ra0_deg = 20.0
    dec0_deg = -30.0
    length_sec = 3600.0 * 8
    freqs = [50]

    # Define telescope settings dictionary
    telescope_settings = get_telescope_settings(
        data_dir,
        ra0_deg,
        dec0_deg,
        length_sec,
        num_time_steps=600,
        swap_xy=True,
    )

    # Copy the base settings dictionary
    current_settings = copy.deepcopy(telescope_settings)

    # Loop over frequencies.
    for freq in freqs:
        # Create the settings file.
        open(settings_path, "w").close()
        settings = oskar.SettingsTree(app1, settings_path)
        settings.from_dict(current_settings)

        # Update output root path and frequency.
        settings["observation/start_frequency_hz"] = 1e6 * freq
        settings["interferometer/ms_filename"] = "testdata.ms"
        settings["sky/oskar_sky_model/file"] = os.path.join(src_dir, "sky.osm")

        # Run the app with the settings file.
        subprocess.call([app1, settings_path])


if __name__ == "__main__":
    main()

    # Add beam info to testdata.ms
    data_dir = os.environ["DATA_DIR"]
    subprocess.check_call(
        [
            "add_beaminfo.py",
            "testdata.ms",
            os.path.join(data_dir, "skalowmini-coef.tm"),
        ]
    )
