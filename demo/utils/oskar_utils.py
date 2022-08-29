#!/usr/bin/env python3
# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
from astropy.time import Time, TimeDelta
from astropy.io import fits
import oskar
import os

"""
Common functions used in the demo project.
"""


def check_tolerance(tolerance, A, B):
    """
    Checks that the maximum (absolute) difference between two matrices are no bigger than the tolerance.

    Parameters
    ----------
    tolerance : float
        (absolute) tolerance
    A : np.ndarray
        Matrix A
    B : np.ndarray
        Matrix B
    """
    if tolerance:
        np.testing.assert_allclose(
            A,
            B,
            atol=tolerance,
            equal_nan=True,
            err_msg=f"Difference between OSKAR and EveryBeam spherical wave model exceeds tolerance of {tolerance}",
        )


def create_settings(app, settings_path, current_settings):
    """
    Create settings file

    app :
        OSKAR application
    settings_path : str
        Path to settings file
    current_settings : dict
        Dictionary with settings
    """
    open(settings_path, "w").close()
    settings = oskar.SettingsTree(app, settings_path)
    settings.from_dict(current_settings)
    return settings


def get_start_time(ra0_deg, length_sec, optimal_time="2000-01-01 12:00:00"):
    """
    Returns optimal start time for field RA and observation length

    Parameters
    ----------
    ra0_deg : float
        Right Ascension (in degrees)
    length_sec : float
        Observation in length
    optimal_time : str, optional
        Optimal time, by default "2000-01-01 12:00:00"

    Returns
    -------
    float
        Start time
    """
    t = Time(optimal_time, scale="utc", location=("116.764d", "0d"))
    dt_hours = (24.0 - t.sidereal_time("apparent").hour) / 1.0027379
    dt_hours += ra0_deg / 15.0
    start = t + TimeDelta(dt_hours * 3600.0 - length_sec / 2.0, format="sec")
    return start.value


def get_telescope_settings(
    data_dir, ra0_deg, dec0_deg, length_sec, num_time_steps, swap_xy
):
    """
    Generate a dict with telescope settings.

    Parameters
    ----------
    data_dir : str
        Path to data directory
    ra0_deg : float
        Right ascension of observation
    dec0_deg : float
        Declination of observation
    length_sec : float
        Observation duration in seconds.
    num_time_steps : int
        Number of time steps
    swap_xy : bool
        Swap x- and y-axis in element coordinate system?

    Returns
    -------
    dict
        Dictionary with telescope settings.
    """
    # Define base settings dictionary.
    return {
        "observation": {
            "phase_centre_ra_deg": ra0_deg,
            "phase_centre_dec_deg": dec0_deg,
            "start_time_utc": get_start_time(
                ra0_deg, length_sec, "2000-01-01 12:00:00"
            ),
            "length": length_sec,
            "num_time_steps": num_time_steps,
        },
        "telescope": {
            "normalise_beams_at_phase_centre": False,
            "aperture_array/element_pattern/normalise": False,
        },
        "telescope/input_directory": os.path.join(
            data_dir, "skalowmini-coef.tm"
        ),
        "telescope/aperture_array/element_pattern/enable_numerical": True,
        "telescope/aperture_array/element_pattern/swap_xy": swap_xy,
        "telescope/aperture_array/array_pattern/enable": True,
        "telescope/aperture_array/array_pattern/normalise": True,
    }


def generate_oskar_csv(
    basefunction_idx, freqs, output_dir, element_id, em_idx=0
):
    """
    Write OSKAR results to csv file.
    """

    l_max = int(np.sqrt(basefunction_idx + 1))

    # Parse all freqs
    for freq in freqs:
        A = np.zeros((l_max, 2 * l_max + 1, 2, 2), dtype=np.complex128)
        data = np.zeros((l_max * (l_max + 2), 4), dtype=np.complex128)
        data[basefunction_idx, em_idx] = 1.0

        i = 0
        for l in range(l_max):
            n = (l + 1) * 2 + 1
            A[l, :n, :] = data[i : i + n, :].reshape(n, 2, 2)
            i += n

        for i, pol in enumerate(("x", "y")):
            for j, em in enumerate(("e", "m")):
                for reim, s in zip(("re", "im"), (1.0, 1.0j)):
                    filename = f"element_pattern_spherical_wave_{pol}_t{em}_{reim}_{element_id}_{freq}.txt"
                    alpha = np.real(A[:, :, i, j] / s)

                    np.savetxt(
                        os.path.join(output_dir, filename),
                        alpha,
                        delimiter=", ",
                    )


def read_oskar_beams(root_name):
    """
    Read oskar beam and store into 4-d numpy array.

    Parameters
    ----------
    root_name : str
        Root name of beam file.

    Returns
    -------
    np.array
        (N, N, 2, 2) array, where N number of image pixels
    """
    A = None
    for i, pol1 in enumerate(["X", "Y"]):
        for j, pol2 in enumerate(["X", "Y"]):
            filename_amp = f"{root_name}_MHz_S0000_TIME_SEP_CHAN_SEP_AMP_{pol1}{pol2}.fits"
            filename_phase = f"{root_name}_MHz_S0000_TIME_SEP_CHAN_SEP_PHASE_{pol1}{pol2}.fits"
            d = (
                fits.getdata(filename_amp)
                * np.exp(1j * fits.getdata(filename_phase))
            )[0, 0, ...].T
            if A is None:
                N = d.shape[-1]
                A = np.zeros((N, N, 2, 2), dtype=np.complex128)

            # One axis in the OSKAR fits files runs in opposite direction compared to the
            # numpy arrays made by make_element_response_image.cpp
            # we flip that axis here to make them the same
            A[..., i, j] = d[::-1, :]

        N = A.shape[0]

    for i in range(N):
        x = 2.0 * i / (N - 1) - 1.0
        for j in range(N):
            y = 2.0 * j / (N - 1) - 1.0
            phi = np.arctan2(y, x)

            e_theta = np.array([[np.cos(phi)], [np.sin(phi)]])
            e_phi = np.array([[-np.sin(phi)], [np.cos(phi)]])

            # This matrix applies the transformation defined in
            # https://github.com/OxfordSKA/OSKAR/blob/master/oskar/convert/define_convert_ludwig3_to_theta_phi_components.h
            T = np.concatenate((e_theta, e_phi), axis=1)

            # Transformation is applied to the right side of A
            # so the rows a of A are tranformed
            A[i, j, ...] = np.dot(A[i, j, ...], T)
    return A
