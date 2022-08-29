#!/usr/bin/env python3
# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
    Script converts LOBEs simulation data to a coefficient file, by fitting spherical harmonics to
    the simulated data. In order to run this script, please make sure
    that the everybeam shared library is on your LD_LIBRARY_PATH and on your PYTHONPATH, See below for
    example export commands:
    - export LD_LIBRARY_PATH=~/opt/everybeam/lib:$LD_LIBRARY_PATH
    - export PYTHONPATH=/home~/opt/everybeam/lib/python3.6/site-packages:$PYTHONPATH

    Script can handle two ways of providing the input data:

    1. station coordinates and simulated results for all frequencies are provided in one HDF5 file
    2. simulated results are provided in a separate MATLAB file per frequency. Station coordinates
    are provided in a separate MATLAB file or retrieved using the lofarantpos tool.
"""

import numpy as np
from scipy.constants import speed_of_light
from scipy.io import loadmat
from everybeam.lobes import F4far_new
import h5py
import os, glob


def read_lofar_antenna_positions(station):
    """
    Get antenna positions in local p,q coordinate system [m]
    within LOFAR station using lofarantpos tool.

    Please note that this function requires lofarantpos to be installed.
    https://pypi.org/project/lofarantpos/

    Parameters
    ----------
    station : str
        LOFAR station name

    Returns
    -------
    np.2darray
        Numpy array (2, nantenna) of p,q-coordinates [m]
    """

    from lofarantpos.db import LofarAntennaDatabase

    db = LofarAntennaDatabase()
    return db.antenna_pqr(station)[:, :2].T


def read_coordinates_from_file(path):
    """
    Read p,q-coordinates from coordinate file at
    specified input path.

    Parameters
    ----------
    path : str
        Path to directory containing simulated results

    Returns
    -------
    np.2darray
        2D array with coordinates in local east-north system [m]
    """

    coord_file = glob.glob(f"{path}/*_coords.mat")
    assert len(coord_file) == 1
    coords = loadmat(coord_file[0])
    p = coords["LBA"]["P"][0, 0].flatten()
    q = coords["LBA"]["Q"][0, 0].flatten()
    return np.array([p, q])


def read_simulation_files(path):
    """
    Read simulation files, where it is assumed that
    each simulation file contains the results for one
    frequency.

    Parameters
    ----------
    path : str
        Path to directory containing simulated results

    Returns
    -------
    dict
        Dictionary, containing the following fields:
        - theta: zenith angles of simulated results [rad]
        - phi: elevation angle of simulated results [rad]
        - freq, np.1darray of frequencies [MHz]
        - v_pol1, returned shape is (#antenna elements, 2, #theta, #phi, #freqs)
        - v_pol2, returned shape is (#antenna elements, 2, #theta, #phi, #freqs)
    """

    file_list = []
    for file_path in glob.glob(f"{path}/*mhz.mat"):
        file_list.append(file_path)

    nfreqs = len(file_list)
    freqs = np.zeros(nfreqs)

    # Hard coded for now, could be inferred when reading the files
    # Alternatively, use np.stack?!
    v_pol1 = np.empty((96, 2, 91, 73, nfreqs), dtype=np.complex128)
    v_pol2 = np.empty((96, 2, 91, 73, nfreqs), dtype=np.complex128)

    for i, file_path in enumerate(file_list):
        print(f"Reading file path {file_path}: {i+1}/{nfreqs}")
        simdata = loadmat(file_path)
        freqs[i] = simdata["Freq"][0, 0]

        if i == 0:
            # Convert to radians
            theta = simdata["Theta"] / 180.0 * np.pi
            phi = simdata["Phi"] / 180.0 * np.pi
        else:
            np.testing.assert_allclose(theta, simdata["Theta"] / 180.0 * np.pi)
            np.testing.assert_allclose(phi, simdata["Phi"] / 180.0 * np.pi)

        v_pol1[..., i] = simdata["Vdiff_pol1"]
        v_pol2[..., i] = simdata["Vdiff_pol2"]

    if np.any(freqs[:-1] > freqs[1:]):
        # Then we have to reorder frequencies - and simulated data -
        # ascendingly
        idcs = np.argsort(freqs)
        freqs = freqs[idcs]
        v_pol1 = v_pol1[..., idcs]
        v_pol2 = v_pol2[..., idcs]

    theta = theta.flatten()
    phi = phi.flatten()

    # Return as dictionary
    # Note: the phi axis runs from [0, 2pi> (not including 2 pi) by cutting off
    # the last element. This avoids 0 and 2pi to be weighted twice.
    return {
        "theta": theta,
        "phi": phi[:-1],
        "v_pol1": v_pol1[..., :-1, :],
        "v_pol2": v_pol2[..., :-1, :],
        "freq": freqs,
    }


def read_h5_file(sim_path, sim_file, freq_select=None, element_select=None):
    """
    Read h5 file with LOBEs simulated data. Assumes the following fields to be present:
        - "X" : east coordinate of antenna element [m]
        - "Y" : north coordinate of antenna element [m]
        - "Freq": frequencies [MHz]
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
            - freq, np.1darray of frequencies [MHz]
            - v_pol1, returned shape is (#antenna elements, 2, #theta, #phi, #freqs)
            - v_pol2, returned shape is (#antenna elements, 2, #theta, #phi, #freqs)
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
            v_pol1_ref = f["Vdiff_pol1"]
            v_pol2_ref = f["Vdiff_pol2"]
        else:
            v_pol1_ref = f["Vdiff_pol1"][freq_select, ..., element_select]
            v_pol2_ref = f["Vdiff_pol2"][freq_select, ..., element_select]

        # Compose complex valued array and transpose
        v_pol1 = (v_pol1_ref["real"] + 1j * v_pol1_ref["imag"]).transpose(
            (4, 3, 2, 1, 0)
        )
        v_pol2 = (v_pol2_ref["real"] + 1j * v_pol2_ref["imag"]).transpose(
            (4, 3, 2, 1, 0)
        )

    # Note: the phi axis runs from [0, 2pi> (not including 2 pi) by cutting off
    # the last element. This avoids 0 and 2pi to be weighted twice.
    return {
        "coords": coords,
        "theta": theta,
        "phi": phi[:-1],
        "freq": freq,
        "v_pol1": v_pol1[..., :-1, :],
        "v_pol2": v_pol2[..., :-1, :],
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


def apply_phasor(v_pol1, v_pol2, freqs, positions, x, y):
    """
    Apply phasor to simulated electromagnetic field.

    Parameters
    ----------
    v_pol1 : np.ndarray
        Differential voltage X-polarization,
        (#elements, 2, #theta, #phi, #freqs)
    v_pol2 : np.ndarray
        Differential voltage X-polarization,
        (#elements, 2, #theta, #phi, #freqs)
    freqs : np.1darray
        Frequencies [Hz]
    positions : np.2darray
        Local (p,q)-coordinates of elements [m]
    x : np.1darray
        Cartesian x-coordinate of (theta, phi) system
    y : np.1darray
        Cartesion y-coordinate of (theta, phi) system
    """

    K = 2.0 * np.pi * freqs / speed_of_light
    dirs = np.array([x, y])

    distance_diff = np.tensordot(positions, dirs, ((0,), (0,)))
    phase_diff = np.multiply.outer(distance_diff, K)
    phasor = np.exp(-1j * phase_diff)

    v_pol1[:] *= phasor[:, np.newaxis, ...]
    v_pol2[:] *= phasor[:, np.newaxis, ...]


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
        Number/order of spherical harmonics
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
        1d array of zenith angle coordinates [rad]
    phi : np.1darray
        1d array of elevation angle coordinates [rad]
    vmin : float, optional
        Minimum value of color range
    vmax : float, optional
        Maximum value of color range
    """
    import matplotlib

    matplotlib.use("tkagg")
    from matplotlib import pyplot as plt

    theta_mesh, phi_mesh = np.meshgrid(theta, phi)
    r = np.sin(theta_mesh)
    xm = r * np.cos(phi_mesh)
    ym = r * np.sin(phi_mesh)

    plt.pcolor(xm, ym, E, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.show()


def fit_and_write_lobes_coeffs(
    theta,
    phi,
    freq,
    positions,
    v_pol1,
    v_pol2,
    nmax,
    outfile,
    apply_weighting=False,
):
    """
    Fit and write coefficients for simulated lobes data

    Parameters
    ----------
    theta : np.1darray
        1d array of zenith angle coordinates [rad]
    phi : np.1darray
        1d array of elevation angle coordinates [rad]
    freq : np.1darray
        Array with frequencies [Hz]
    positions : np.2darray
        Element coordinates in local east-north system [m]
    v_pol1 : np.5darray
        Simulated voltage pattern for X-polarisation, is expected to have
        shape (#antenna elements, 2, #theta, #phi, #freqs)
    v_pol2 : np.5darray
        Simulated voltage pattern for X-polarisation, is expected to have
        shape (#antenna elements, 2, #theta, #phi, #freqs)
    nmax : int
        Number/order of spherical harmonics
    outfile : str
        Path to output h5 file
    apply_weighting : bool, optional
        Apply weighting with the sine of the zenith angle, i.e. sin(theta), by default False

    Returns
    -------
    np.2darray, np.ndarray
        Return basis function expansion and coefficients to ease further postprocessing
    """

    # Frequencies are in MHz, convert to Hz
    sintheta = np.sin(theta.reshape(-1, theta.size))
    weights = None if apply_weighting is None else sintheta.flatten()

    cosphi = np.cos(phi.reshape(-1, phi.size))
    sinphi = np.sin(phi.reshape(-1, phi.size))

    x = sintheta.T * cosphi
    y = sintheta.T * sinphi

    # NOTE: Always check that freqs are in Hz
    apply_phasor(v_pol1, v_pol2, freq, positions, x, y)

    # X/Y-polarization
    E1, E1w = extract_electromagnetic_field(v_pol1, weights=weights)
    E2, E2w = extract_electromagnetic_field(v_pol2, weights=weights)

    F, F_w, nms = compute_polyomial_expansion(
        theta, phi, nmax, weights=weights
    )
    print("Inverting matrix")
    Finv = np.linalg.pinv(F_w) if F_w is not None else np.linalg.pinv(F)

    print("Compute coefficints as Finv \dot E")
    p1 = np.tensordot(Finv, E1w if E1w is not None else E1, axes=([1], [2]))
    p2 = np.tensordot(Finv, E2w if E2w is not None else E2, axes=([1], [2]))

    # Stack into 4D array
    p = np.stack((p1, p2), axis=-1)
    write_h5file(outfile, p.T, nms, freq)
    return F, p.T


def main():
    """
    There are two routes for fitting coefficients to the
    LOBEs simulations:
    - Route 1: station coordinates and simulated results for all frequencies are provided in one HDF5 file
    - Route 2: simulated results are provided in a separate HDF5 file per frequency. Also, station coordinates
    are provided in a separate HDF5 file.

    Examples for both routes - with paths and filenames on my local system - are provided below.
    """

    # Order of Legendre polynomial to be fitted. Adjust to suit your needs
    nmax = 21
    approach = "separate_hdf5_file"

    def single_hdf5_file():
        """
        Coordinates and simulated results for all frequencies in one HDF5 file
        """

        # Adjust to local paths
        station = "SE607"
        simpath = "/var/scratch/maljaars/data"
        simfile = f"{station}_patterns.mat"
        outfile = f"LOBES_{station}LBA.h5"
        outfile = os.path.join(".", outfile)
        freq_select = None
        element_select = None

        h5data = read_h5_file(
            simpath,
            simfile,
            freq_select=freq_select,
            element_select=element_select,
        )
        # Please note that freqs have to be converted from MHz to Hz!
        F, pT, = fit_and_write_lobes_coeffs(
            h5data["theta"],
            h5data["phi"],
            h5data["freq"] * 1e6,
            h5data["coords"],
            h5data["v_pol1"],
            h5data["v_pol2"],
            nmax,
            outfile,
            apply_weighting=False,
        )
        return simdata["theta"], simdata["phi"], F, pT

    def separate_hdf5_files():
        """
        Simulated results for each frequency stored in separate file. Station coordinates
        stored in separate HDF5 file.
        """

        # Adjust to correct local paths
        station = "CS007"
        path = f"/var/scratch/maljaars/data/LOBES/ftp.astron.nl/outgoing/arts/{station}"
        outpath = "/var/scratch/maljaars/data/LOBES/core-station-coeffs"
        outfile = f"LOBES_{station}LBA.h5"
        outfile = os.path.join(outpath, outfile)

        # Read positions from file
        # positions = read_coordinates_from_file(path)
        # Use lofarantpos to read coordinates
        positions = read_lofar_antenna_positions(f"{station}LBA")
        simdata = read_simulation_files(path)

        F, pT, = fit_and_write_lobes_coeffs(
            simdata["theta"],
            simdata["phi"],
            simdata["freq"],
            positions,
            simdata["v_pol1"],
            simdata["v_pol2"],
            nmax,
            outfile,
            apply_weighting=False,
        )
        return simdata["theta"], simdata["phi"], F, pT

    theta, phi, F, pT = (
        single_hdf5_file()
        if approach == "single_hdf5_file"
        else separate_hdf5_files()
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
