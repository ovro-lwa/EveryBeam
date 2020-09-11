#!/usr/bin/env python3
"""
Generates all-sky zenith-centred beam patterns for SKALA-4 and EDA-2 antennas.
"""

import copy
import glob
import os
import re
import subprocess

from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
# pylint: disable=wrong-import-position
import matplotlib.pyplot as plt
import numpy
from utils import create_settings, get_start_time


def main(npixels):
    """Main function."""
    # Name of the application to run, and a settings file for it.
    app1 = 'oskar_sim_interferometer'
    app2 = 'oskar_sim_beam_pattern'
    settings_path = '_temp_settings.ini'

    data_dir = os.environ["DATA_DIR"]

    # Define some basic observation parameters.
    ra0_deg = 0.0
    dec0_deg = -27.0
    length_sec = 60 # 3600.0*6

    # Define base settings dictionary.
    common_settings = {
        'observation': {
            'phase_centre_ra_deg': ra0_deg,
            'phase_centre_dec_deg': dec0_deg,
            #'pointing_file': 'station_pointing.txt',
            'start_time_utc': get_start_time(ra0_deg, length_sec, '2000-01-01 12:00:00'),
            'length': length_sec,
            'num_time_steps': 60

        },
        'telescope': {
            'normalise_beams_at_phase_centre': False,
            'aperture_array/element_pattern/normalise': False
        }
    }

    # Define frequencies of interest (in MHz).
    freqs = [50]
    #freqs = [70]

    # Define telescope models to use, and associated overrides for them.
    telescopes = {
        #'stationresponse': {
        'basefunctions': {
            'telescope/input_directory': os.path.join(data_dir, 'skalowmini-coef.tm'),
            'telescope/aperture_array/element_pattern/enable_numerical': True,
            'telescope/aperture_array/element_pattern/swap_xy': False,
            'telescope/aperture_array/array_pattern/enable': True,
            'telescope/aperture_array/array_pattern/normalise': True
        },
    }

    # Loop over telescope models.
    for tel, tel_params in telescopes.items():

        # Copy the base settings dictionary.
        current_settings = copy.deepcopy(common_settings)

        # Update current settings with telescope model parameters.
        current_settings.update(tel_params)

        # Loop over frequencies.
        for freq in freqs:

            # Create the settings file.
            settings = create_settings(app1, settings_path, current_settings)

            # Update output root path and frequency.
            tel_root = re.sub(r'[^\w]', '', tel)  # Strip symbols from tel.
            root_path = tel_root + ('_%03d_MHz' % freq)
            settings['observation/start_frequency_hz'] = 1e6 * freq
            settings['interferometer/ms_filename'] = 'skalowmini-coef.MS'

            # Run the app with the settings file.
            subprocess.call([app1, settings_path])

            # Create the settings file.
            settings = create_settings(app2, settings_path, current_settings)

            # Update output root path and frequency.
            tel_root = re.sub(r'[^\w]', '', tel)  # Strip symbols from tel.
            root_path = tel_root + ('_%03d_MHz' % freq)
            settings['beam_pattern/root_path'] = root_path
            settings['observation/start_frequency_hz'] = 1e6 * freq
            settings['beam_pattern/coordinate_frame'] = 'Horizon'
            settings['beam_pattern/beam_image/size'] = npixels
            settings['beam_pattern/station_outputs/fits_image/amp'] = True
            settings['beam_pattern/station_outputs/fits_image/phase'] = True
            settings['beam_pattern/station_outputs/fits_image/auto_power_real'] = False

            # Run the app with the settings file.
            subprocess.call([app2, settings_path])


if __name__ == '__main__':
    main(256)
