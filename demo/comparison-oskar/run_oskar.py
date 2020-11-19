#!/usr/bin/env python3
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Generates all-sky zenith-centred beam patterns for SKALA-4 and EDA-2 antennas.
"""

import copy
import glob
import os.path
import re
import subprocess

from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
# pylint: disable=wrong-import-position
import matplotlib.pyplot as plt
import numpy
from utils import create_settings, get_start_time


def plot_panel(ax, image, title, cmap):
    """Plots a single panel."""
    im = ax.imshow(numpy.squeeze(image), cmap=cmap)
    plt.colorbar(im, format='%.2e')
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    ax.set_title(title)
    ax.axis('equal')


def make_plots(title, glob_pattern, out_basename):
    """Generates a plot with four panels."""
    # Load FITS images matching the glob pattern.
    files = glob.glob(glob_pattern)
    files.sort()  # Must be sorted!
    images = []
    for file in files:
        images.append(fits.getdata(file))

    # Set titles and colour maps to use.
    cmap = ''
    titles = []
    if '_AUTO_POWER' in glob_pattern:
        # cmap = 'plasma'
        cmap = 'Blues_r'
        titles = ['Stokes I', 'Stokes Q', 'Stokes U', 'Stokes V']
    elif '_AMP' in glob_pattern:
        # cmap = 'jet'
        cmap = 'CMRmap'
        titles = ['XX', 'XY', 'YX', 'YY']

    # Sub-plots, one for each Stokes parameter.
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(221, frameon=False)
    plot_panel(ax, images[0], titles[0], cmap)
    ax = fig.add_subplot(222, frameon=False)
    plot_panel(ax, images[1], titles[1], cmap)
    ax = fig.add_subplot(223, frameon=False)
    plot_panel(ax, images[2], titles[2], cmap)
    ax = fig.add_subplot(224, frameon=False)
    plot_panel(ax, images[3], titles[3], cmap)

    # Add main title.
    title = title.replace('_', ' ')  # Replace underscores with spaces.
    fig.suptitle(title)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)

    # Save and close.
    plt.savefig('%s.png' % (out_basename))
    plt.close('all')


def main(npixels):
    """Main function."""
    # Name of the application to run, and a settings file for it.
    app = 'oskar_sim_beam_pattern'
    settings_path = '_temp_settings.ini'

    # Define some basic observation parameters.
    ra0_deg = 0.0
    dec0_deg = -27.0
    length_sec = 1.0

    # Define base settings dictionary.
    common_settings = {
        'observation': {
            'phase_centre_ra_deg': ra0_deg,
            'phase_centre_dec_deg': dec0_deg,
            #'pointing_file': 'station_pointing.txt',
            'start_time_utc': get_start_time(ra0_deg, length_sec, '2000-01-01 00:00:00'),
            'length': length_sec
        },
        'telescope': {
            'normalise_beams_at_phase_centre': False,
            'aperture_array/element_pattern/normalise': False
        },
        'beam_pattern': {
            'coordinate_frame': 'Horizon',
            'beam_image/size': npixels,
            'station_outputs/fits_image/amp': True,
            'station_outputs/fits_image/phase': True,
            'station_outputs/fits_image/auto_power_real': False
        }
    }

    # Define frequencies of interest (in MHz).
    freqs = [50]
    #freqs = [70]

    # Define telescope models to use, and associated overrides for them.
    telescopes = {
        'basefunctions': {
            'telescope/input_directory': 'telescope.tm',
            'telescope/aperture_array/element_pattern/enable_numerical': True,
            'telescope/aperture_array/element_pattern/swap_xy': False,
            'telescope/aperture_array/array_pattern/enable': False
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
            settings = create_settings(app, settings_path, current_settings)

            # Update output root path and frequency.
            tel_root = re.sub(r'[^\w]', '', tel)  # Strip symbols from tel.
            root_path = tel_root + ('_%03d_MHz' % freq)
            settings['beam_pattern/root_path'] = root_path
            settings['observation/start_frequency_hz'] = 1e6 * freq

            # Run the app with the settings file.
            subprocess.call([app, settings_path])

            # Make plots.
            #title = tel + ' @ ' + str(freq) + ' MHz'
            #make_plots(title=title,
                       #glob_pattern=root_path+'*_AMP*',
                       #out_basename=root_path+'_amp')
            #make_plots(title=title,
                       #glob_pattern=root_path+'*_AUTO_POWER*',
                       #out_basename=root_path+'_stokes')


if __name__ == '__main__':
    main(32)
