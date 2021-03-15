#!/usr/bin/env python
# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Script converts LOBEs simulation data to a coefficient file, by projecting the
# simulated data onto spherical harmonics. In order to run this script, please make sure
# that the everybeam shared lybrary is on your LD_LIBRARY_PATH and on your PYTHONPATH, See below for
# example export commands:
# - export LD_LIBRARY_PATH=~/opt/everybeam/lib:$LD_LIBRARY_PATH
# - export PYTHONPATH=/home~/opt/everybeam/lib/python3.6/site-packages:$PYTHONPATH

import numpy as np
from scipy.constants import speed_of_light
from everybeam.lobes import F4far_new
import h5py
import os


def read_h5_file(sim_path, sim_file, freq_select=None, element_select=None):
    """
    Read h5 file with LOBEs simulated data. Assumes the following fields to be present:
        - "X" : east coordinate of antenna element [m]
        - "Y" : north coordinate of antenna element [m]
        - "Theta" : zenith angles for which the embedded element patterns are calculated [deg]
        - "Phi" : elevation angles for which the embeeded element patterns are calculated [deg]
        - "Vdiff_pol1" : differential voltage for polarization 1, assumed shape (#freqs, #phi, #theta, 2, #antenna elements)
        - "Vdiff_pol2" : differential voltage for polarization 2, assumed shape (#freqs, #phi, #theta, 2, #antenna elements)

    Parameters
    ----------
    sim_path : str
        Path containing the simulated data
    sim_file : str
        File name containing the simulated data

    Returns
    -------
    dict
        Dictionary containing numpy arrays of the simulated data. Contains the following keys:
            - coords
            - theta, zenith angle [rad]
            - phi, elevation angle [rad]
            - freq
            - V_pol1, returned shape is (#antenna elements, 2, #theta, #phi, #freqs)
            - V_pol1, returned shape is (#antenna elements, 2, #theta, #phi, #freqs)
    """
    with h5py.File(os.path.join(sim_path, sim_file), "r") as f:
        X = f["X"] if element_select is None else f["X"][:, element_select]
        Y = f["Y"] if element_select is None else f["Y"][:, element_select]
        coords = np.concatenate((X, Y), axis=0)

        theta = np.deg2rad(np.array(f["Theta"]).flatten())
        phi = np.deg2rad(np.array(f["Phi"]).flatten())

        freq = (
            np.array(f["Freq"]).flatten()
            if freq_select is None
            else np.array(f["Freq"]).flatten()[freq_select]
        )

        if freq_select is None and element_select is None:
            V_pol1_ref = f["Vdiff_pol1"]
            V_pol2_ref = f["Vdiff_pol2"]
        else:
            V_pol1_ref = f["Vdiff_pol1"][freq_select, ..., element_select]
            V_pol2_ref = f["Vdiff_pol2"][freq_select, ..., element_select]

        # Compose complex valued array and transpose
        V_pol1 = (V_pol1_ref["real"] + 1j * V_pol1_ref["imag"]).transpose(
            (4, 3, 2, 1, 0)
        )
        V_pol2 = (V_pol2_ref["real"] + 1j * V_pol2_ref["imag"]).transpose(
            (4, 3, 2, 1, 0)
        )

    # Note: the phi axis runs from [0, 2pi> (not including 2 pi) by cutting off
    # the last element. This avoids 0 and 2pi to be weighted twice.
    return {
        "coords": coords,
        "theta": theta,
        "phi": phi[:-1],
        "freq": freq,
        "V_pol1": V_pol1[..., :-1, :],
        "V_pol2": V_pol2[..., :-1, :],
    }


def extract_electromagnetic_field(V_pol, weights=None):
    """
    Extract the simulated electromagnetic field.

    Parameters
    ----------
    V_pol : np.ndarray
        Differential voltage, shape expected to be
        (#elements, 2, #theta, #phi, #freqs)
    weights : np.1darray
        1darray with (zenith angle) weights with size equal to #theta, by default None

    Returns
    -------
    np.ndarray, np.ndarray/None
        Unweighted electromagnetic field, Weighted electromagnetic field (None if weights is None).
        Shape equals (#element, #freqs, #phi * #theta)
    """
    nr_elements = V_pol.shape[0]
    nr_frequencies = V_pol.shape[-1]

    Etheta = V_pol[:, 0, ...].transpose((0, 3, 2, 1))
    Ephi = V_pol[:, 1, ...].transpose((0, 3, 2, 1))

    if weights is not None:
        Etheta_w = Etheta * weights
        Ephi_w = Ephi * weights

        Etheta_w = Etheta_w.reshape((nr_elements, nr_frequencies, -1))
        Ephi_w = Ephi_w.reshape((nr_elements, nr_frequencies, -1))

    Etheta = Etheta.reshape((nr_elements, nr_frequencies, -1))
    Ephi = Ephi.reshape((nr_elements, nr_frequencies, -1))

    if weights is not None:
        return (
            np.concatenate((Etheta, Ephi), axis=2),
            np.concatenate((Etheta_w, Ephi_w), axis=2),
        )
    else:
        return np.concatenate((Etheta, Ephi), axis=2), None


def compute_polyomial_expansion(theta, phi, nmax, weights=None):
    """
    Compute the (Legendre) polynomial expansion

    Parameters
    ----------
    theta : np.1darray
        1d array of zenith angles [rad]
    phi : np.1darray
        1d array of elevation angles [rad]
    nmax : int
        Number of spherical harmonics
    weights : np.1darray, optional
        1d of zenith angle weights, should have same size as theta array.
        Defaults to None, in which case no weight

    Returns
    -------
    np.2darray, np.2darray, np.ndarray
        Basis matrix, Weighted basis matrix, nms matrix. If weights is None,
        the weighted basis matrix will be None
    """

    Theta, Phi = np.meshgrid(theta, phi)
    print(f"Theta size {Theta.size}")
    F = np.zeros((2 * Theta.size, 2 * nmax * (nmax + 2)), dtype=np.complex128)

    jn = 0
    nms = []
    for n in range(1, nmax + 1):
        for m in range(-n, n + 1):
            for s in [1, 2]:
                q2, q3 = F4far_new(s, m, n, Theta.flatten(), Phi.flatten())
                F[:, jn] = np.concatenate((q2, q3))
                nms.append((n, m, s))
                jn += 1

    if weights is not None:
        weights = np.tile(weights, 2 * phi.size)
        F_w = F * weights.reshape(-1, 1)
        return F, F_w, nms
    else:
        return F, None, nms


def write_h5file(outfile, coeffs, nms, freqs):
    """
    Write results to h5 file

    Parameters
    ----------
    outfile : str
        Output file name (or path)
    coeffs : np.ndarray
        5d numpy array, expected dimensions (2, #freqs, #elements, #coefficients)
    nms : np.ndarray
    freqs : np.1darray
        1D array with frequencies
    outdir : str, optional
        Output directory, by default "."
    """
    with h5py.File(outfile, "w") as f:
        f.create_dataset("coefficients", data=coeffs)
        f.create_dataset("nms", data=nms)
        f.create_dataset("frequencies", data=freqs.astype(np.float64))


def make_plot(E, theta, phi, vmin=None, vmax=None):
    """
    Make plot of field E on x, y plane

    Parameters
    ----------
    E : np.ndarray
        Electromagnetic field
    theta : np.1darray
        1d array of zenith angle coordinates (in rad)
    phi : np.1darray
        1d array of elevation angle
    vmin : float, optional
        Minimum value of color range
    vmax : float, optional
        Maximum value of color range
    """
    import matplotlib

    matplotlib.use("tkagg")
    from matplotlib import pyplot as plt

    Theta, Phi = np.meshgrid(theta, phi)
    r = np.sin(Theta)
    xm = r * np.cos(Phi)
    ym = r * np.sin(Phi)

    plt.pcolor(xm, ym, E, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.show()


def fit_and_write_lobes_coeffs(
    simpath,
    simfile,
    outfile,
    nmax,
    freq_select=None,
    element_select=None,
    apply_weighting=False,
):
    """
    Fit and write coefficients on simulated lobes data

    Parameters
    ----------
    simpath : str
        Path where simulated data file can be found
    simfile : str
        File name of simulated data file
    outfile : str
        Output file
    nmax: int
        Number of spherical harmonics
    freq_select : slice, optional
        Slicer object, can be used to select only certain frequency bands, by default None in which case
        all frequencies are fitted. Use with care, validity of slicer is not checked!
    element_select : slice, optional
        Slicer object, can be used to select only certain elements, by default None in which case
        all frequencies are fitted. Use with care, validity of slicer is not checked!
    apply_weighting : bool, optional
        Apply weighting with the sine of the zenith angle, i.e. sin(theta), by default False

    Returns
    -------
    np.1darray, np.1darray, np.2darray, np.ndarray
        Return theta, phi, basis function expansion and coefficients to ease further postprocessing
    """

    h5data = read_h5_file(
        simpath, simfile, freq_select=freq_select, element_select=element_select
    )

    theta = h5data["theta"]
    phi = h5data["phi"]
    freq = h5data["freq"]
    positions = h5data["coords"]
    V_pol1 = h5data["V_pol1"]
    V_pol2 = h5data["V_pol2"]

    # Frequencies are in MHz, convert to Hz
    K = 2.0 * np.pi * freq * 1e6 / speed_of_light
    sintheta = np.sin(theta.reshape(-1, theta.size))
    weights = None if apply_weighting is None else sintheta.flatten()

    cosphi = np.cos(phi.reshape(-1, phi.size))
    sinphi = np.sin(phi.reshape(-1, phi.size))

    x = sintheta.T * cosphi
    y = sintheta.T * sinphi

    dirs = np.array([x, y])
    distance_diff = np.tensordot(positions, dirs, ((0,), (0,)))
    phase_diff = np.multiply.outer(distance_diff, K)

    phasor = np.exp(-1j * phase_diff)
    V_pol1 *= phasor[:, np.newaxis, ...]
    V_pol2 *= phasor[:, np.newaxis, ...]

    # X/Y-polarization
    E1, E1w = extract_electromagnetic_field(V_pol1, weights=weights)
    E2, E2w = extract_electromagnetic_field(V_pol2, weights=weights)

    F, F_w, nms = compute_polyomial_expansion(theta, phi, nmax, weights=weights)
    print("Inverting matrix")
    Finv = np.linalg.pinv(F_w) if F_w is not None else np.linalg.pinv(F)

    print("Compute coefficints as Finv \dot E")
    p1 = np.tensordot(Finv, E1w if E1w is not None else E1, axes=([1], [2]))
    p2 = np.tensordot(Finv, E2w if E2w is not None else E2, axes=([1], [2]))

    # Stack into 4D array
    p = np.stack((p1, p2), axis=-1)
    write_h5file(outfile, p.T, nms, freq)
    return theta, phi, F, p.T


def main():
    # Adjust these to suit your needs
    simpath = "/var/scratch/maljaars/data"
    simfile = "SE607_patterns.mat"
    outfile = "LOBES_SE607LBA.h5"
    outfile = os.path.join(".", outfile)
    nmax = 21
    freq_select = None
    element_select = None
    theta, phi, F, pT, = fit_and_write_lobes_coeffs(
        simpath,
        simfile,
        outfile,
        nmax,
        freq_select=freq_select,
        element_select=element_select,
        apply_weighting=False,
    )

    # Plot modulus of fitted response for XX, freq 0, and element 0
    make_plot(
        np.abs(
            np.dot(F[0 : theta.size * phi.size, :], pT[0, 0, 0, :]).reshape(
                phi.size, theta.size
            )
        ),
        theta,
        phi,
    )


if __name__ == "__main__":
    main()
