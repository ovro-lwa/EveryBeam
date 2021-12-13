// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_CIRCULARSYMMETRIC_ATCACOEFFICIENTS_H_
#define EVERYBEAM_CIRCULARSYMMETRIC_ATCACOEFFICIENTS_H_

#include <array>
#include <map>
#include <string>

#include "coefficients.h"

namespace everybeam {
namespace circularsymmetric {
class ATCACoefficients final : public Coefficients {
 public:
  aocommon::UVector<double> GetFrequencies(
      [[maybe_unused]] double frequency) const override {
    constexpr size_t n_bands = 7;
    constexpr double first_band_frequency = 1332e6;
    constexpr double frequency_step = 256e6;
    aocommon::UVector<double> frequencies(n_bands);
    for (size_t i = 0; i != n_bands; ++i) {
      frequencies[i] = first_band_frequency + i * frequency_step;
    }
    return frequencies;
  }
  aocommon::UVector<double> GetCoefficients(
      [[maybe_unused]] double frequency) const override {
    return aocommon::UVector<double>(std::begin(coefficients_),
                                     std::end(coefficients_));
  }
  double MaxRadiusInArcMin() const override { return 53.0; }
  double ReferenceFrequency() const override { return 1e9; }
  bool AreInverted() const override { return true; }

 private:
  // coef x nfreq
  static constexpr std::array<double, 35> coefficients_{
      1.0, 1.06274e-03, 1.32342e-06, -8.72013e-10, 1.08020e-12,
      1.0, 9.80817e-04, 1.17898e-06, -7.83160e-10, 8.66199e-13,
      1.0, 9.53553e-04, 9.33233e-07, -4.26759e-10, 5.63667e-13,
      1.0, 9.78268e-04, 6.63231e-07, 4.18235e-11,  2.62297e-13,
      1.0, 1.02424e-03, 6.12726e-07, 2.25733e-10,  2.04834e-13,
      1.0, 1.05818e-03, 5.37473e-07, 4.22386e-10,  1.17530e-13,
      1.0, 1.10650e-03, 5.11574e-07, 5.89732e-10,  8.13628e-14};
};
}  // namespace circularsymmetric
}  // namespace everybeam
#endif  // EVERYBEAM_CIRCULARSYMMETRIC_ATCACOEFFICIENTS_H_
