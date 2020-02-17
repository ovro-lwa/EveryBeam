#!/usr/bin/python3

import argparse
import os
import scipy.io as sio
from os.path import dirname, join as pjoin
import numpy as np
import h5py
from tqdm import tqdm

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_dir", type=str, help="directory containing coefficients")
parser.add_argument("output_file", type=str, help="name of the HDF5 output file")
parser.add_argument("-lmax", type=int, help="select the number of coefficients in one dimension")
args = parser.parse_args()

# Read the input directory
# This directory should contain subdirectories,
# with names corresponding to frequency.
input_dir = args.input_dir
print("Reading coefficients from: " , input_dir)
dirs = os.listdir(input_dir)
dirs = list(filter(lambda x: x.isdigit(), dirs))
freqs = sorted(map(int, dirs))
print("Parsing frequencies %d-%d MHz" % (int(freqs[0]), int(freqs[-1])))

# Select the upper-left rectangle of size l_max * l_max
l_max = args.lmax

# Create HDF5 file
output_file = args.output_file
print("Creating HDF5 file: ", output_file)
hf = h5py.File(output_file, 'a')

# Parse all freqs
for freq in tqdm(freqs):
    #print("Parsing frequency: %d MHz" % freq)

    # Get path to matlab data files
    patha = pjoin(input_dir, str(freq), 'alphas_pola.mat')
    pathb = pjoin(input_dir, str(freq), 'alphas_polb.mat')

    # Read the matlab data files
    matrixa = sio.loadmat(patha)
    matrixb = sio.loadmat(pathb)

    # Convert to numpy arrays
    alpha_te_a = np.asarray(matrixa['alpha_te'])
    alpha_tm_a = np.asarray(matrixa['alpha_tm'])
    alpha_te_b = np.asarray(matrixb['alpha_te'])
    alpha_tm_b = np.asarray(matrixb['alpha_tm'])

    # Make sure that all four sub-matrices are shaped identically
    assert(alpha_te_a.shape == alpha_tm_a.shape == alpha_te_b.shape == alpha_tm_b.shape)
    shape = alpha_te_a.shape
    nr_elements = shape[0]
    nr_coeffs_y = shape[1]
    nr_coeffs_x = shape[2]

    # Select coeffs
    if (l_max and
        l_max < nr_coeffs_y and
        l_max < nr_coeffs_y):
        alpha_te_a = alpha_te_a[:,:l_max,:l_max]
        alpha_tm_a = alpha_tm_a[:,:l_max,:l_max]
        alpha_te_b = alpha_te_b[:,:l_max,:l_max]
        alpha_tm_b = alpha_tm_b[:,:l_max,:l_max]
        shape = alpha_te_a.shape

    # Put coeffs in one matrix
    data = np.array([alpha_te_a, alpha_tm_a, alpha_te_b, alpha_tm_b])
    # [tetm*pol,ant,l,m]

    # Swap antenna and polarization*tetm axis
    # [ant,pol*tetm,tetm,l,m]
    data = np.swapaxes(data, 0, 1)

    # Split tetm and polarization axis
    # [ant,pol,tetm,tetm,l,m]
    s = data.shape
    assert(s[1] == 4)
    s = (s[0],) + (2,2) + s[2:]
    data = data.reshape(s)

    # Sanity check
    assert(np.allclose(data[:, 0, 0, :, :], alpha_te_a))
    assert(np.allclose(data[:, 1, 0, :, :], alpha_te_b))
    assert(np.allclose(data[:, 1, 1, :, :], alpha_tm_b))
    assert(np.allclose(data[:, 0, 1, :, :], alpha_tm_a))

    # Add group for current frequency
    hf.create_dataset(str(freq), s, data=data)

# Close HDF5 file
hf.close()

# Done
print("Done!")
