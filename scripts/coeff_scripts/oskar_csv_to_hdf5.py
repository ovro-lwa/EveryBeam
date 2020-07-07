#!/usr/bin/env python3

import argparse
import os
import scipy.io as sio
import os.path
import numpy as np
import h5py
from tqdm import tqdm
from glob import glob
import csv

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

# Only read element_id 0, because EveryBeam's implementation of the OSKAR model
# currently supports only a single model for all elements
element_id = 0
nr_elements = 1

files = glob(os.path.join(input_dir, 'element_pattern_spherical_wave_?_t?_??_{}_*.txt').format(element_id))

files = [os.path.basename(f) for f in files]

freqs = sorted(set([int(f[:-4].split('_')[-1]) for f in files]))


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

    A = None

    for i, pol in enumerate(('x', 'y')):
        for j, em in enumerate(('e', 'm')):
            for reim, s in zip(('re', 'im'), (1.0, 1.0j)) :
                filename = 'element_pattern_spherical_wave_{}_t{}_{}_{}_{}.txt'.format(pol, em, reim, element_id, freq)
                alpha = np.loadtxt(os.path.join(input_dir,filename), delimiter=',', ndmin=2)
                l_max = alpha.shape[0]
                assert(alpha.shape[1] == 2*l_max + 1)
                if A is None:
                    A = np.zeros( (l_max,2*l_max+1, 2,2), dtype=np.complex128)
                else:
                    assert(A.shape == (l_max,2*l_max+1, 2,2))
                A[:,:,i,j] += alpha*s

    # Create the output matrix
    data = np.zeros((nr_elements, l_max * (l_max+2), 4), dtype=np.complex128)

    # Fill the output matrix
    # The inner dimension are set as follows:
    # (x_te_re, x_te_im), (x_tm_re, x_tm_im),
    # (y_te_re, y_te_im), (y_tm_re, y_tm_im)

    i = 0
    for l in range(l_max):
        n = (l+1) * 2 + 1
        data[0,i:i+n,:] = A[l,:n,:].reshape(n,4)
        i += n

    # Add group for current frequency
    hf.create_dataset(str(freq), data.shape, data=data)

# Close HDF5 file
hf.close()

# Done
print("Done!")
