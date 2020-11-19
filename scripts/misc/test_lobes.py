#!/usr/bin/env python3
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np

import lobes

a = np.linspace(0, 0.5, 10)

lobes.scale(a, 2.0)

print(lobes.legendre(1,-11,a))
