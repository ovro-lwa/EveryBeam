// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_OSKAR_H_
#define EVERYBEAM_OSKAR_H_

#include "oskar_vector_types.h"
#include "oskar_helper.h"

void oskar_evaluate_dipole_pattern_double(const int num_points,
                                          const double* theta,
                                          const double* phi,
                                          const double freq_hz,
                                          const double dipole_length_m,
                                          std::complex<double>* pattern);

void oskar_evaluate_spherical_wave_sum_double(double theta, double phi_x,
                                              double phi_y, int l_max,
                                              const std::complex<double>* alpha,
                                              std::complex<double>* pattern);

#endif