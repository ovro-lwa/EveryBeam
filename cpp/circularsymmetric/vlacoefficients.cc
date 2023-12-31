// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "vlacoefficients.h"

#include <cmath>
#include <vector>

using everybeam::circularsymmetric::VLACoefficients;

std::array<double, 5> VLACoefficients::GetCoefficients(
    const std::string& band_name, double freq) {
  char band = '?';

  const size_t sharp = band_name.find('#');
  if (sharp != std::string::npos) {
    if (sharp > 5 && band_name.substr(0, 5) == "EVLA_") band = band_name[5];
  }
  if (band == '?') {
    band = DetermineFeed(freq);
  }
  LimitFreqForBand(band, freq);

  auto coeffmap = GetCoefficients();

  const std::array<double, 5>* coeff;
  const double freqMHz = freq * 1e-6;
  std::map<int, std::array<double, 5>>::iterator low, prev;
  low = coeffmap.lower_bound(freqMHz);
  if (low == coeffmap.end()) {
    --low;
    coeff = &low->second;
  } else if (low == coeffmap.begin()) {
    coeff = &low->second;
  } else {
    prev = low;
    --prev;
    if (std::fabs(freqMHz - prev->first) < std::fabs(low->first - freqMHz)) {
      coeff = &prev->second;
    } else {
      coeff = &low->second;
    }
  }
  return *coeff;
}

char VLACoefficients::DetermineFeed(double freq, double freq_center) {
  if ((freq_center > 224e6 && freq_center < 480e6) ||
      (freq > 224e6 && freq < 480e6)) {
    return 'P';
  }
  if ((freq_center > 900e6 && freq_center < 2003.0e6) ||
      (freq > 900e6 && freq < 2003e6)) {
    return 'L';
  }
  if ((freq_center > 1990e6 && freq_center < 4001.0e6) ||
      (freq > 1990e6 && freq < 4001e6)) {
    return 'S';
  }
  if ((freq_center > 3990e6 && freq_center < 8001.0e6) ||
      (freq > 3990e6 && freq < 8001e6)) {
    return 'C';
  }
  if ((freq_center > 7990e6 && freq_center < 12001.0e6) ||
      (freq > 7990e6 && freq < 12001e6)) {
    return 'X';
  }
  if ((freq_center > 12000e6 && freq_center < 18000.0e6) ||
      (freq > 12000e6 && freq < 18000e6)) {
    return 'U';
  }
  if ((freq_center > 19000e6 && freq_center < 26000.0e6) ||
      (freq > 19000e6 && freq < 26000e6)) {
    return 'K';
  }
  if ((freq_center > 28000e6 && freq_center < 38000.0e6) ||
      (freq > 28000e6 && freq < 38000e6)) {
    return 'A';
  }
  if ((freq_center > 41000e6 && freq_center < 50000.0e6) ||
      (freq > 41000e6 && freq < 50000e6)) {
    return 'Q';
  }
  return '?';
}

// From PBMath1DEVLA::limitFreqForBand
void VLACoefficients::LimitFreqForBand(char band, double& freq) {
  if (band == 'P') {
    if (freq <= 224e6) freq = 232e6;
    if (freq >= 480e6) freq = 470e6;
  } else if (band == 'L') {
    if (freq <= 900e6) freq = 1040e6;
    if (freq >= 2000e6) freq = 2000e6;
  } else if (band == 'S') {
    if (freq < 2052e6) freq = 2052e6;
    if (freq >= 3948e6) freq = 3948e6;
  } else if (band == 'C') {
    if (freq < 4052e6) freq = 4052e6;
    if (freq >= 7948e6) freq = 7948e6;
  } else if (band == 'X') {
    if (freq < 8052e6) freq = 8052e6;
    if (freq >= 11948e6) freq = 11948e6;
  } else if (band == 'U') {
    if (freq < 12052e6) freq = 12052e6;
    if (freq >= 17948e6) freq = 17948e6;
  } else if (band == 'K') {
    if (freq < 19052e6) freq = 19052e6;
    if (freq >= 25948e6) freq = 25948e6;
  } else if (band == 'A') {
    if (freq < 28052e6) freq = 28052e6;
    if (freq >= 38048e6) freq = 38048e6;
  } else if (band == 'Q') {
    if (freq < 41052e6) freq = 41052e6;
    if (freq >= 43948e6) freq = 43948e6;
  }
}

std::map<char, double> VLACoefficients::GetFeedConf() {
  std::map<char, double> feed_conf;
  feed_conf['L'] =
      (-185.9) * M_PI / 180.0;  // squint orientation, rads, North of +AZ axis
  feed_conf['S'] = (-11.61) * M_PI / 180.0;
  feed_conf['C'] = (-104.8) * M_PI / 180.0;
  feed_conf['X'] = (-113.7) * M_PI / 180.0;
  feed_conf['U'] = (42.4) * M_PI / 180.0;
  feed_conf['K'] = (64.4) * M_PI / 180.0;
  feed_conf['A'] = (106.9) * M_PI / 180.0;
  feed_conf['Q'] = (85.5) * M_PI / 180.0;
  return feed_conf;
}

std::map<int, std::array<double, 5>> VLACoefficients::GetCoefficients() {
  // This comes from PBMath1DEVLA::init()
  // Also see https://library.nrao.edu/public/memos/evla/EVLAM_195.pdf

  /*
   * For reference, these are the correponding frequencies in MHz:
   *
   * {232., 246., 281., 296., 312., 328., 344., 357.,
   *  382., 392., 403., 421., 458., 470., 1040, 1104,
   *  1168, 1232, 1296, 1360, 1424, 1488, 1552, 1680,
   *  1744, 1808, 1872, 1936, 2000};
   */

  std::map<int, std::array<double, 5>> coeffmap;
  ////P
  coeffmap[232] = {1.0, -1.137e-3, 5.19e-7, -1.04e-10, 0.71e-14};
  coeffmap[246] = {1.0, -1.130e-3, 5.04e-7, -1.02e-10, 0.77e-14};
  coeffmap[281] = {1.0, -1.106e-3, 5.11e-7, -1.10e-10, 0.91e-14};
  coeffmap[296] = {1.0, -1.125e-3, 5.27e-7, -1.14e-10, 0.96e-14};
  coeffmap[312] = {1.0, -1.030e-3, 4.44e-7, -0.89e-10, 0.68e-14};
  coeffmap[328] = {1.0, -0.980e-3, 4.25e-7, -0.87e-10, 0.69e-14};
  coeffmap[344] = {1.0, -0.974e-3, 4.09e-7, -0.76e-10, 0.53e-14};
  coeffmap[357] = {1.0, -0.996e-3, 4.23e-7, -0.79e-10, 0.51e-14};
  coeffmap[382] = {1.0, -1.002e-3, 4.39e-7, -0.88e-10, 0.64e-14};
  coeffmap[392] = {1.0, -1.067e-3, 5.13e-7, -1.12e-10, 0.90e-14};
  coeffmap[403] = {1.0, -1.057e-3, 4.90e-7, -1.06e-10, 0.87e-14};
  coeffmap[421] = {1.0, -1.154e-3, 5.85e-7, -1.33e-10, 1.08e-14};
  coeffmap[458] = {1.0, -0.993e-3, 4.67e-7, -1.04e-10, 0.88e-14};
  coeffmap[470] = {1.0, -1.010e-3, 4.85e-7, -1.07e-10, 0.86e-14};
  /////////L
  coeffmap[1040] = {1.000, -1.529e-3, 8.69e-7, -1.88e-10};
  coeffmap[1104] = {1.000, -1.486e-3, 8.15e-7, -1.68e-10};
  coeffmap[1168] = {1.000, -1.439e-3, 7.53e-7, -1.45e-10};
  coeffmap[1232] = {1.000, -1.450e-3, 7.87e-7, -1.63e-10};
  coeffmap[1296] = {1.000, -1.428e-3, 7.62e-7, -1.54e-10};
  coeffmap[1360] = {1.000, -1.449e-3, 8.02e-7, -1.74e-10};
  coeffmap[1424] = {1.000, -1.462e-3, 8.23e-7, -1.83e-10};
  coeffmap[1488] = {1.000, -1.455e-3, 7.92e-7, -1.63e-10};
  coeffmap[1552] = {1.000, -1.435e-3, 7.54e-7, -1.49e-10};
  coeffmap[1680] = {1.000, -1.443e-3, 7.74e-7, -1.57e-10};
  coeffmap[1744] = {1.000, -1.462e-3, 8.02e-7, -1.69e-10};
  coeffmap[1808] = {1.000, -1.488e-3, 8.38e-7, -1.83e-10};
  coeffmap[1872] = {1.000, -1.486e-3, 8.26e-7, -1.75e-10};
  coeffmap[1936] = {1.000, -1.459e-3, 7.93e-7, -1.62e-10};
  coeffmap[2000] = {1.000, -1.508e-3, 8.31e-7, -1.68e-10};
  ////////S
  coeffmap[2052] = {1.000, -1.429e-3, 7.52e-7, -1.47e-10};
  coeffmap[2180] = {1.000, -1.389e-3, 7.06e-7, -1.33e-10};
  coeffmap[2436] = {1.000, -1.377e-3, 6.90e-7, -1.27e-10};
  coeffmap[2564] = {1.000, -1.381e-3, 6.92e-7, -1.26e-10};
  coeffmap[2692] = {1.000, -1.402e-3, 7.23e-7, -1.40e-10};
  coeffmap[2820] = {1.000, -1.433e-3, 7.62e-7, -1.54e-10};
  coeffmap[2948] = {1.000, -1.433e-3, 7.46e-7, -1.42e-10};
  coeffmap[3052] = {1.000, -1.467e-3, 8.05e-7, -1.70e-10};
  coeffmap[3180] = {1.000, -1.497e-3, 8.38e-7, -1.80e-10};
  coeffmap[3308] = {1.000, -1.504e-3, 8.37e-7, -1.77e-10};
  coeffmap[3436] = {1.000, -1.521e-3, 8.63e-7, -1.88e-10};
  coeffmap[3564] = {1.000, -1.505e-3, 8.37e-7, -1.75e-10};
  coeffmap[3692] = {1.000, -1.521e-3, 8.51e-7, -1.79e-10};
  coeffmap[3820] = {1.000, -1.534e-3, 8.57e-7, -1.77e-10};
  coeffmap[3948] = {1.000, -1.516e-3, 8.30e-7, -1.66e-10};
  /// C
  coeffmap[4052] = {1.000, -1.406e-3, 7.41e-7, -1.48e-10};
  coeffmap[4180] = {1.000, -1.385e-3, 7.09e-7, -1.36e-10};
  coeffmap[4308] = {1.000, -1.380e-3, 7.08e-7, -1.37e-10};
  coeffmap[4436] = {1.000, -1.362e-3, 6.95e-7, -1.35e-10};
  coeffmap[4564] = {1.000, -1.365e-3, 6.92e-7, -1.31e-10};
  coeffmap[4692] = {1.000, -1.339e-3, 6.56e-7, -1.17e-10};
  coeffmap[4820] = {1.000, -1.371e-3, 7.06e-7, -1.40e-10};
  coeffmap[4948] = {1.000, -1.358e-3, 6.91e-7, -1.34e-10};
  coeffmap[5052] = {1.000, -1.360e-3, 6.91e-7, -1.33e-10};
  coeffmap[5180] = {1.000, -1.353e-3, 6.74e-7, -1.25e-10};
  coeffmap[5308] = {1.000, -1.359e-3, 6.82e-7, -1.27e-10};
  coeffmap[5436] = {1.000, -1.380e-3, 7.05e-7, -1.37e-10};
  coeffmap[5564] = {1.000, -1.376e-3, 6.99e-7, -1.31e-10};
  coeffmap[5692] = {1.000, -1.405e-3, 7.39e-7, -1.47e-10};
  coeffmap[5820] = {1.000, -1.394e-3, 7.29e-7, -1.45e-10};
  coeffmap[5948] = {1.000, -1.428e-3, 7.57e-7, -1.57e-10};
  coeffmap[6052] = {1.000, -1.445e-3, 7.68e-7, -1.50e-10};
  coeffmap[6148] = {1.000, -1.422e-3, 7.38e-7, -1.38e-10};
  coeffmap[6308] = {1.000, -1.463e-3, 7.94e-7, -1.62e-10};
  coeffmap[6436] = {1.000, -1.478e-3, 8.22e-7, -1.74e-10};
  coeffmap[6564] = {1.000, -1.473e-3, 8.00e-7, -1.62e-10};
  coeffmap[6692] = {1.000, -1.455e-3, 7.76e-7, -1.53e-10};
  coeffmap[6820] = {1.000, -1.487e-3, 8.22e-7, -1.72e-10};
  coeffmap[6948] = {1.000, -1.472e-3, 8.05e-7, -1.67e-10};
  coeffmap[7052] = {1.000, -1.470e-3, 8.01e-7, -1.64e-10};
  coeffmap[7180] = {1.000, -1.503e-3, 8.50e-7, -1.84e-10};
  coeffmap[7308] = {1.000, -1.482e-3, 8.19e-7, -1.72e-10};
  coeffmap[7436] = {1.000, -1.498e-3, 8.22e-7, -1.66e-10};
  coeffmap[7564] = {1.000, -1.490e-3, 8.18e-7, -1.66e-10};
  coeffmap[7692] = {1.000, -1.481e-3, 7.98e-7, -1.56e-10};
  coeffmap[7820] = {1.000, -1.474e-3, 7.94e-7, -1.57e-10};
  coeffmap[7948] = {1.000, -1.448e-3, 7.69e-7, -1.51e-10};
  //////X
  coeffmap[8052] = {1.000, -1.403e-3, 7.21e-7, -1.37e-10};
  coeffmap[8180] = {1.000, -1.398e-3, 7.10e-7, -1.32e-10};
  coeffmap[8308] = {1.000, -1.402e-3, 7.16e-7, -1.35e-10};
  coeffmap[8436] = {1.000, -1.400e-3, 7.12e-7, -1.32e-10};
  coeffmap[8564] = {1.000, -1.391e-3, 6.95e-7, -1.25e-10};
  coeffmap[8692] = {1.000, -1.409e-3, 7.34e-7, -1.49e-10};
  coeffmap[8820] = {1.000, -1.410e-3, 7.36e-7, -1.45e-10};
  coeffmap[8948] = {1.000, -1.410e-3, 7.34e-7, -1.43e-10};
  coeffmap[9052] = {1.000, -1.403e-3, 7.20e-7, -1.36e-10};
  coeffmap[9180] = {1.000, -1.396e-3, 7.09e-7, -1.31e-10};
  coeffmap[9308] = {1.000, -1.432e-3, 7.68e-7, -1.55e-10};
  coeffmap[9436] = {1.000, -1.414e-3, 7.43e-7, -1.47e-10};
  coeffmap[9564] = {1.000, -1.416e-3, 7.45e-7, -1.47e-10};
  coeffmap[9692] = {1.000, -1.406e-3, 7.26e-7, -1.39e-10};
  coeffmap[9820] = {1.000, -1.412e-3, 7.36e-7, -1.43e-10};
  coeffmap[9948] = {1.000, -1.409e-3, 7.29e-7, -1.39e-10};
  coeffmap[10052] = {1.000, -1.421e-3, 7.46e-7, -1.45e-10};
  coeffmap[10180] = {1.000, -1.409e-3, 7.25e-7, -1.36e-10};
  coeffmap[10308] = {1.000, -1.402e-3, 7.13e-7, -1.31e-10};
  coeffmap[10436] = {1.000, -1.399e-3, 7.09e-7, -1.29e-10};
  coeffmap[10564] = {1.000, -1.413e-3, 7.37e-7, -1.43e-10};
  coeffmap[10692] = {1.000, -1.412e-3, 7.34e-7, -1.41e-10};
  coeffmap[10820] = {1.000, -1.401e-3, 7.12e-7, -1.31e-10};
  coeffmap[10948] = {1.000, -1.401e-3, 7.12e-7, -1.31e-10};
  coeffmap[10052] = {1.000, -1.401e-3, 7.12e-7, -1.31e-10};
  coeffmap[11180] = {1.000, -1.394e-3, 6.99e-7, -1.24e-10};
  coeffmap[11308] = {1.000, -1.394e-3, 7.01e-7, -1.26e-10};
  coeffmap[11436] = {1.000, -1.391e-3, 6.94e-7, -1.22e-10};
  coeffmap[11564] = {1.000, -1.389e-3, 6.92e-7, -1.22e-10};
  coeffmap[11692] = {1.000, -1.386e-3, 6.80e-7, -1.15e-10};
  coeffmap[11820] = {1.000, -1.391e-3, 6.88e-7, -1.19e-10};
  coeffmap[11948] = {1.000, -1.399e-3, 6.97e-7, -1.22e-10};
  /// U
  coeffmap[12052] = {1.000, -1.399e-3, 7.17e-7, -1.34e-10};
  coeffmap[12180] = {1.000, -1.392e-3, 7.07e-7, -1.31e-10};
  coeffmap[12308] = {1.000, -1.393e-3, 7.19e-7, -1.38e-10};
  coeffmap[12436] = {1.000, -1.393e-3, 7.20e-7, -1.40e-10};
  coeffmap[12564] = {1.000, -1.395e-3, 7.19e-7, -1.38e-10};
  coeffmap[12692] = {1.000, -1.397e-3, 7.20e-7, -1.37e-10};
  coeffmap[12820] = {1.000, -1.388e-3, 7.06e-7, -1.32e-10};
  coeffmap[12948] = {1.000, -1.397e-3, 7.18e-7, -1.36e-10};
  coeffmap[13052] = {1.000, -1.400e-3, 7.27e-7, -1.40e-10};
  coeffmap[13180] = {1.000, -1.406e-3, 7.44e-7, -1.50e-10};
  coeffmap[13308] = {1.000, -1.403e-3, 7.37e-7, -1.47e-10};
  coeffmap[13436] = {1.000, -1.392e-3, 7.08e-7, -1.31e-10};
  coeffmap[13564] = {1.000, -1.384e-3, 6.94e-7, -1.24e-10};
  coeffmap[13692] = {1.000, -1.382e-3, 6.95e-7, -1.25e-10};
  coeffmap[13820] = {1.000, -1.376e-3, 6.88e-7, -1.24e-10};
  coeffmap[13948] = {1.000, -1.384e-3, 6.98e-7, -1.28e-10};
  coeffmap[14052] = {1.000, -1.400e-3, 7.36e-7, -1.48e-10};
  coeffmap[14180] = {1.000, -1.397e-3, 7.29e-7, -1.45e-10};
  coeffmap[14308] = {1.000, -1.399e-3, 7.32e-7, -1.45e-10};
  coeffmap[14436] = {1.000, -1.396e-3, 7.25e-7, -1.42e-10};
  coeffmap[14564] = {1.000, -1.393e-3, 7.20e-7, -1.39e-10};
  coeffmap[14692] = {1.000, -1.384e-3, 7.03e-7, -1.31e-10};
  coeffmap[14820] = {1.000, -1.388e-3, 7.06e-7, -1.32e-10};
  coeffmap[14948] = {1.000, -1.393e-3, 7.16e-7, -1.37e-10};
  coeffmap[15052] = {1.000, -1.402e-3, 7.38e-7, -1.48e-10};
  coeffmap[15180] = {1.000, -1.407e-3, 7.47e-7, -1.53e-10};
  coeffmap[15308] = {1.000, -1.406e-3, 7.41e-7, -1.48e-10};
  coeffmap[15436] = {1.000, -1.399e-3, 7.31e-7, -1.44e-10};
  coeffmap[15564] = {1.000, -1.397e-3, 7.28e-7, -1.43e-10};
  coeffmap[15692] = {1.000, -1.401e-3, 7.35e-7, -1.46e-10};
  coeffmap[15820] = {1.000, -1.402e-3, 7.34e-7, -1.45e-10};
  coeffmap[15948] = {1.000, -1.399e-3, 7.30e-7, -1.44e-10};
  coeffmap[16052] = {1.000, -1.419e-3, 7.59e-7, -1.54e-10};
  coeffmap[16180] = {1.000, -1.419e-3, 7.59e-7, -1.52e-10};
  coeffmap[16308] = {1.000, -1.412e-3, 7.40e-7, -1.44e-10};
  coeffmap[16436] = {1.000, -1.407e-3, 7.32e-7, -1.40e-10};
  coeffmap[16564] = {1.000, -1.408e-3, 7.32e-7, -1.41e-10};
  coeffmap[16692] = {1.000, -1.410e-3, 7.34e-7, -1.40e-10};
  coeffmap[16820] = {1.000, -1.407e-3, 7.27e-7, -1.38e-10};
  coeffmap[16948] = {1.000, -1.423e-3, 7.63e-7, -1.55e-10};
  coeffmap[17052] = {1.000, -1.437e-3, 7.87e-7, -1.66e-10};
  coeffmap[17180] = {1.000, -1.438e-3, 7.84e-7, -1.64e-10};
  coeffmap[17308] = {1.000, -1.445e-3, 7.98e-7, -1.71e-10};
  coeffmap[17436] = {1.000, -1.452e-3, 8.10e-7, -1.77e-10};
  coeffmap[17564] = {1.000, -1.458e-3, 8.13e-7, -1.70e-10};
  coeffmap[17692] = {1.000, -1.456e-3, 8.06e-7, -1.72e-10};
  coeffmap[17820] = {1.000, -1.453e-3, 8.00e-7, -1.68e-10};
  coeffmap[17948] = {1.000, -1.452e-3, 7.99e-7, -1.69e-10};
  /////K
  coeffmap[19052] = {1.000, -1.419e-3, 7.56e-7, -1.53e-10};
  coeffmap[19180] = {1.000, -1.426e-3, 7.70e-7, -1.59e-10};
  coeffmap[19308] = {1.000, -1.433e-3, 7.82e-7, -1.64e-10};
  coeffmap[19436] = {1.000, -1.429e-3, 7.73e-7, -1.60e-10};
  coeffmap[19564] = {1.000, -1.427e-3, 7.70e-7, -1.59e-10};
  coeffmap[19692] = {1.000, -1.425e-3, 7.65e-7, -1.56e-10};
  coeffmap[19820] = {1.000, -1.430e-3, 7.76e-7, -1.62e-10};
  coeffmap[19948] = {1.000, -1.434e-3, 7.81e-7, -1.63e-10};
  coeffmap[21052] = {1.000, -1.448e-3, 8.05e-7, -1.73e-10};
  coeffmap[21180] = {1.000, -1.436e-3, 7.84e-7, -1.63e-10};
  coeffmap[21308] = {1.000, -1.441e-3, 7.94e-7, -1.68e-10};
  coeffmap[21436] = {1.000, -1.439e-3, 7.89e-7, -1.66e-10};
  coeffmap[21564] = {1.000, -1.442e-3, 7.96e-7, -1.69e-10};
  coeffmap[21692] = {1.000, -1.435e-3, 7.81e-7, -1.61e-10};
  coeffmap[21820] = {1.000, -1.442e-3, 7.92e-7, -1.66e-10};
  coeffmap[21948] = {1.000, -1.439e-3, 7.82e-7, -1.61e-10};
  coeffmap[23052] = {1.000, -1.401e-3, 7.21e-7, -1.37e-10};
  coeffmap[23180] = {1.000, -1.408e-3, 7.31e-7, -1.41e-10};
  coeffmap[23308] = {1.000, -1.407e-3, 7.28e-7, -1.39e-10};
  coeffmap[23436] = {1.000, -1.407e-3, 7.31e-7, -1.41e-10};
  coeffmap[23564] = {1.000, -1.419e-3, 7.47e-7, -1.47e-10};
  coeffmap[23692] = {1.000, -1.395e-3, 7.10e-7, -1.33e-10};
  coeffmap[23820] = {1.000, -1.413e-3, 7.36e-7, -1.42e-10};
  coeffmap[23948] = {1.000, -1.402e-3, 7.21e-7, -1.36e-10};
  coeffmap[25052] = {1.000, -1.402e-3, 7.17e-7, -1.31e-10};
  coeffmap[25180] = {1.000, -1.432e-3, 7.73e-7, -1.58e-10};
  coeffmap[25308] = {1.000, -1.407e-3, 7.22e-7, -1.33e-10};
  coeffmap[25436] = {1.000, -1.417e-3, 7.43e-7, -1.45e-10};
  coeffmap[25564] = {1.000, -1.422e-3, 7.52e-7, -1.48e-10};
  coeffmap[25692] = {1.000, -1.427e-3, 7.59e-7, -1.52e-10};
  coeffmap[25820] = {1.000, -1.416e-3, 7.42e-7, -1.44e-10};
  coeffmap[25948] = {1.000, -1.422e-3, 7.46e-7, -1.45e-10};
  /// A
  coeffmap[28052] = {1.000, -1.444e-3, 7.61e-7, -1.44e-10};
  coeffmap[28180] = {1.000, -1.439e-3, 7.54e-7, -1.42e-10};
  coeffmap[28308] = {1.000, -1.457e-3, 7.87e-7, -1.58e-10};
  coeffmap[28436] = {1.000, -1.457e-3, 7.90e-7, -1.60e-10};
  coeffmap[28564] = {1.000, -1.455e-3, 7.87e-7, -1.59e-10};
  coeffmap[28692] = {1.000, -1.458e-3, 7.88e-7, -1.58e-10};
  coeffmap[28820] = {1.000, -1.453e-3, 7.81e-7, -1.56e-10};
  coeffmap[28948] = {1.000, -1.460e-3, 7.98e-7, -1.64e-10};
  coeffmap[31052] = {1.000, -1.415e-3, 7.44e-7, -1.44e-10};
  coeffmap[31180] = {1.000, -1.408e-3, 7.26e-7, -1.37e-10};
  coeffmap[31308] = {1.000, -1.413e-3, 7.28e-7, -1.36e-10};
  coeffmap[31436] = {1.000, -1.394e-3, 7.07e-7, -1.30e-10};
  coeffmap[31564] = {1.000, -1.404e-3, 7.23e-7, -1.37e-10};
  coeffmap[31692] = {1.000, -1.427e-3, 7.48e-7, -1.44e-10};
  coeffmap[31820] = {1.000, -1.418e-3, 7.48e-7, -1.48e-10};
  coeffmap[31948] = {1.000, -1.413e-3, 7.37e-7, -1.42e-10};
  coeffmap[34052] = {1.000, -1.42e-3, 7.28e-7, -1.34e-10};
  coeffmap[34180] = {1.000, -1.46e-3, 7.77e-7, -1.53e-10};
  coeffmap[34308] = {1.000, -1.42e-3, 7.41e-7, -1.42e-10};
  coeffmap[34436] = {1.000, -1.42e-3, 7.36e-7, -1.39e-10};
  coeffmap[34564] = {1.000, -1.46e-3, 7.76e-7, -1.52e-10};
  coeffmap[34692] = {1.000, -1.42e-3, 7.34e-7, -1.38e-10};
  coeffmap[34820] = {1.000, -1.42e-3, 7.34e-7, -1.39e-10};
  coeffmap[34948] = {1.000, -1.45e-3, 7.68e-7, -1.49e-10};
  coeffmap[37152] = {1.000, -1.42e-3, 7.47e-7, -1.44e-10};
  coeffmap[37280] = {1.000, -1.41e-3, 7.35e-7, -1.40e-10};
  coeffmap[37408] = {1.000, -1.45e-3, 7.65e-7, -1.46e-10};
  coeffmap[37536] = {1.000, -1.41e-3, 7.13e-7, -1.29e-10};
  coeffmap[37664] = {1.000, -1.41e-3, 7.30e-7, -1.38e-10};
  coeffmap[37792] = {1.000, -1.45e-3, 7.75e-7, -1.50e-10};
  coeffmap[37820] = {1.000, -1.45e-3, 7.68e-7, -1.49e-10};
  coeffmap[38048] = {1.000, -1.41e-3, 7.38e-7, -1.43e-10};
  // Q
  coeffmap[41052] = {1.000, -1.453e-3, 7.69e-7, -1.47e-10};
  coeffmap[41180] = {1.000, -1.479e-3, 8.03e-7, -1.61e-10};
  coeffmap[41308] = {1.000, -1.475e-3, 7.97e-7, -1.58e-10};
  coeffmap[41436] = {1.000, -1.451e-3, 7.73e-7, -1.51e-10};
  coeffmap[41564] = {1.000, -1.450e-3, 7.71e-7, -1.51e-10};
  coeffmap[41692] = {1.000, -1.465e-3, 7.79e-7, -1.49e-10};
  coeffmap[41820] = {1.000, -1.460e-3, 7.73e-7, -1.47e-10};
  coeffmap[41948] = {1.000, -1.434e-3, 7.47e-7, -1.40e-10};
  coeffmap[43052] = {1.000, -1.428e-3, 7.40e-7, -1.38e-10};
  coeffmap[43180] = {1.000, -1.418e-3, 7.29e-7, -1.34e-10};
  coeffmap[43308] = {1.000, -1.433e-3, 7.49e-7, -1.43e-10};
  coeffmap[43436] = {1.000, -1.438e-3, 7.55e-7, -1.45e-10};
  coeffmap[43564] = {1.000, -1.419e-3, 7.36e-7, -1.40e-10};
  coeffmap[43692] = {1.000, -1.397e-3, 7.13e-7, -1.33e-10};
  coeffmap[43820] = {1.000, -1.423e-3, 7.39e-7, -1.40e-10};
  coeffmap[43948] = {1.000, -1.452e-3, 7.68e-7, -1.47e-10};
  return coeffmap;
}
