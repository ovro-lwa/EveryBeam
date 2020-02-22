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
parser.add_argument("-debug", dest="debug", action="store_true")
parser.set_defaults(input_dir=".")
parser.set_defaults(output_file="oskar.h5")
parser.set_defaults(debug=False)
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

# Create HDF5 file
output_file = args.output_file
print("Creating HDF5 file: ", output_file)
hf = h5py.File(output_file, 'w')

# Debugging
debug = args.debug
if (debug):
    freqs = [freqs[0]]
    tqdm = lambda x: x

# Parse all freqs
for freq in tqdm(freqs):
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

    # Get dimensions of input matrix
    # This matrix is lower triangular, e.g.:
    # 32x65 (2080 elements), with 1088 (32x(32+2)) nonzeros
    shape = alpha_te_a.shape
    nr_elements = shape[0]
    src_height = shape[1]
    src_width = shape[2]
    if (debug):
        print("nr_elements: ", nr_elements)
        print("input height: ", src_height)
        print("input width: ", src_width)

    # Determine dimensions of output matrix
    # This matrix is rectangular and contains
    # all nonzero elements of the input matrix
    dst_height = src_height + 2
    dst_width = src_height
    if (debug):
        print("output height: ", dst_height)
        print("output width: ", dst_width)

    # Fill the output matrix
    # The inner dimension are set as follows:
    # (x_te_re, x_te_im), (x_tm_re, x_tm_im),
    # (y_te_re, y_te_im), (y_tm_re, y_tm_im)
    data = np.zeros((nr_elements, dst_height * dst_width, 4), dtype=np.complex128)

    i = 0
    for m in range(0,src_height):
        for l in range((m+1)*2, -1, -1):
            ind0 = l * src_width + m
            data[:,i,0] = alpha_te_a[:,m,l]
            data[:,i,1] = alpha_tm_a[:,m,l]
            data[:,i,2] = alpha_te_b[:,m,l]
            data[:,i,3] = alpha_tm_b[:,m,l]
            i += 1

    # Sanity check
    assert(np.count_nonzero(data[:,:,0]) == np.count_nonzero(alpha_te_a))
    assert(np.count_nonzero(data[:,:,1]) == np.count_nonzero(alpha_tm_a))
    assert(np.count_nonzero(data[:,:,2]) == np.count_nonzero(alpha_te_b))
    assert(np.count_nonzero(data[:,:,3]) == np.count_nonzero(alpha_tm_b))

    # Add group for current frequency
    hf.create_dataset(str(freq), data.shape, data=data)

# Close HDF5 file
hf.close()

# Done
print("Done!")