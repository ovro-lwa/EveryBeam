#!/usr/bin/env python3
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse
import os
import scipy.io as sio
import os.path
import numpy as np
import h5py
from tqdm import tqdm
from glob import glob

output_dir = 'telescope.tm'
element_id = 0


freqs = [50]

def generate_oskar_csv(basefunction_idx, em_idx=0):

    l_max = int(np.sqrt(basefunction_idx+1))

    # Parse all freqs
    for freq in freqs:

        A = np.zeros( (l_max,2*l_max+1, 2,2), dtype=np.complex128)
        data = np.zeros((l_max * (l_max+2), 4), dtype=np.complex128)
        data[basefunction_idx,em_idx] = 1.0

        i = 0
        for l in range(l_max):
            n = (l+1) * 2 + 1
            A[l,:n,:] = data[i:i+n,:].reshape(n,2,2)
            i += n

        for i, pol in enumerate(('x', 'y')):
            for j, em in enumerate(('e', 'm')):
                for reim, s in zip(('re', 'im'), (1.0, 1.0j)) :
                    filename = 'element_pattern_spherical_wave_{}_t{}_{}_{}_{}.txt'.format(pol, em, reim, element_id, freq)
                    alpha = np.real(A[:,:,i,j]/s)

                    np.savetxt(os.path.join(output_dir,filename), alpha, delimiter=', ')

if __name__ == "__main__":
    generate_oskar_csv(0)
