#!/usr/bin/env python3
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Generates all-sky zenith-centred beam patterns for SKALA-4 and EDA-2 antennas.
"""

import copy
import subprocess

from oskar_utils import create_settings, get_telescope_settings
from plot_utils import make_fits_plots


def main(npixels):
    """Main function."""
    # Name of the application to run, and a settings file for it.
    app = "oskar_sim_beam_pattern"
    settings_path = "_temp_settings.ini"

    # Define some basic observation parameters.
    ra0_deg = 0.0
    dec0_deg = -27.0
    length_sec = 1.0
    data_dir = "."  # Not relevant
    num_steps = 1  # Not relevant

    # Define frequencies of interest (in MHz).
    freqs = [50]

    # Define settings dictionary.
    common_settings = get_telescope_settings(
        data_dir, ra0_deg, dec0_deg, length_sec, num_steps, False
    )
    beam_pattern = {
        "beam_pattern": {
            "coordinate_frame": "Horizon",
            "beam_image/size": npixels,
            "station_outputs/fits_image/amp": True,
            "station_outputs/fits_image/phase": True,
            "station_outputs/fits_image/auto_power_real": False,
        }
    }
    # Override and extend telescope settings
    telescope = {
        "telescope/input_directory": "telescope.tm",
        "telescope/aperture_array/element_pattern/enable_numerical": True,
        "telescope/aperture_array/element_pattern/swap_xy": False,
        "telescope/aperture_array/array_pattern/enable": False,
    }
    common_settings = {**common_settings, **beam_pattern, **telescope}

    # Copy the base settings dictionary.
    current_settings = copy.deepcopy(common_settings)

    # Loop over frequencies.
    for freq in freqs:
        # Create the settings file.
        settings = create_settings(app, settings_path, current_settings)

        # Update output root path and frequency.
        root_path = "basefunctions" + f"_{int(freq):03}_MHz"
        settings["beam_pattern/root_path"] = root_path
        settings["observation/start_frequency_hz"] = 1e6 * freq

        # Run the app with the settings file.
        subprocess.call([app, settings_path])

        # [OPTIONALLY] Make plots.
        # title = f"basefunctions @ {freq} MHz"
        # make_fits_plots(title=title,
        # glob_pattern=root_path+'*_AMP*',
        # out_basename=root_path+'_amp')
        # make_fits_plots(title=title,
        # glob_pattern=root_path+'*_AUTO_POWER*',
        # out_basename=root_path+'_stokes')


if __name__ == "__main__":
    main(32)
