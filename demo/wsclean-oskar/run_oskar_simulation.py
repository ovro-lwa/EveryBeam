#!/usr/bin/env python3
"""
Generates all-sky zenith-centred beam patterns for SKALA-4 and EDA-2 antennas.
"""

import copy
import glob
import os.path
import re
import subprocess

from astropy.io import fits
from astropy.time import Time, TimeDelta
import matplotlib
matplotlib.use('Agg')
# pylint: disable=wrong-import-position
import matplotlib.pyplot as plt
import numpy
import oskar


def get_start_time(ra0_deg, length_sec):
    """Returns optimal start time for field RA and observation length."""
    t = Time('2000-01-01 12:00:00', scale='utc', location=('116.764d', '0d'))
    dt_hours = (24.0 - t.sidereal_time('apparent').hour) / 1.0027379
    dt_hours += (ra0_deg / 15.0)
    start = t + TimeDelta(dt_hours * 3600.0 - length_sec / 2.0, format='sec')
    return start.value

def main():
    """Main function."""
    # Name of the application to run, and a settings file for it.
    app1 = 'oskar_sim_interferometer'
    settings_path = '_temp_settings.ini'
    data_dir = os.environ["DATA_DIR"]
    src_dir = os.path.dirname(__file__)

    # Define some basic observation parameters.
    ra0_deg = 20.0
    dec0_deg = -30.0
    length_sec = 3600.0*8

    # Define base settings dictionary.
    common_settings = {
        'observation': {
            'phase_centre_ra_deg': ra0_deg,
            'phase_centre_dec_deg': dec0_deg,
            #'pointing_file': 'station_pointing.txt',
            'start_time_utc': get_start_time(ra0_deg, length_sec),
            'length': length_sec,
            'num_time_steps': 600

        },
        'telescope': {
            'normalise_beams_at_phase_centre': False,
            'aperture_array/element_pattern/normalise': False
        },
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
            'telescope/aperture_array/element_pattern/swap_xy': True,
            'telescope/aperture_array/array_pattern/enable': True,
            'telescope/aperture_array/array_pattern/normalise': True
        },
    }

    # Loop over telescope models.
    for tel, tel_params in telescopes.items():

        # Copy the base settings diction.ary.
        current_settings = copy.deepcopy(common_settings)

        # Update current settings with telescope model parameters.
        current_settings.update(tel_params)

        # Loop over frequencies.
        for freq in freqs:

            # Create the settings file.
            open(settings_path, 'w').close()
            settings = oskar.SettingsTree(app1, settings_path)
            settings.from_dict(current_settings)

            # Update output root path and frequency.
            tel_root = re.sub(r'[^\w]', '', tel)  # Strip symbols from tel.
            root_path = tel_root + ('_%03d_MHz' % freq)
            settings['observation/start_frequency_hz'] = 1e6 * freq
            settings['interferometer/ms_filename'] = 'testdata.ms'
            settings['sky/oskar_sky_model/file'] = os.path.join(src_dir, 'sky.osm')

            # Run the app with the settings file.
            subprocess.call([app1, settings_path])

if __name__ == '__main__':
    main()
