// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_CIRCULARSYMMETRIC_GMRTCOEFFICIENTS_H_
#define EVERYBEAM_CIRCULARSYMMETRIC_GMRTCOEFFICIENTS_H_

#include <array>
#include <map>
#include <string>

#include "coefficients.h"

namespace everybeam {
namespace circularsymmetric {
/**
 * References:
 *
 * -
 * http://www.ncra.tifr.res.in/ncra/gmrt/gmrt-users/observing-help/ugmrt-primary-beam-shape
 * - http://www.gmrt.ncra.tifr.res.in/doc/beam-shape-v1-09sep2022.pdf
 * - https://github.com/ruta-k/uGMRTprimarybeam
 * - PBCOR in AIPS and code/synthesis/TransformMachines/PBMath.cc in CASA.
 */
class GMRTCoefficients final : public Coefficients {
 public:
  aocommon::UVector<double> GetFrequencies(
      [[maybe_unused]] double frequency) const override {
    return {187.5e6,  // band 2 (125-250 MHz)
            375e6,    // band 3 (250-500 MHz)
            700e6,    // band 4 (550-850 MHz)
            1250e6};  // band 5 (1050-1450 MHz)
  }
  aocommon::UVector<double> GetCoefficients(
      [[maybe_unused]] double frequency) const override {
    return aocommon::UVector<double>(std::begin(coefficients_),
                                     std::end(coefficients_));
  }
  double MaxRadiusInArcMin() const override {
    constexpr double diameter = 45.0;
    constexpr double frequency = 43.0;
    return 1.1998662 * 25.0 / diameter * frequency;
  }
  double ReferenceFrequency() const override { return 1.5e9; }
  bool AreInverted() const override { return false; }

 private:
  /**
   * These come from the PDF with title
   * "Upgraded GMRT: (Updated) primary-beam shape parameters",
   * http://www.gmrt.ncra.tifr.res.in/doc/beam-shape-v1-09sep2022.pdf
   *
   * 2023-01-30: band 4 en 5 were updated and band 2 was added.
   *
   * This array has dimensions n_coef x n_freq, with n_coef=5 and
   * n_freq=4 (matching with @ref GetFrequencies() ).
   */
  static constexpr std::array<double, 20> coefficients_{
      1.0, -3.089e-03, 39.314e-07, -23.011e-10, 5.037e-13,
      1.0, -2.939e-03, 33.312e-07, -16.659e-10, 3.066e-13,
      1.0, -3.263e-03, 42.618e-07, -25.580e-10, 5.823e-13,
      1.0, -2.614e-03, 27.594e-07, -13.268e-10, 2.395e-13,
  };
};
}  // namespace circularsymmetric
}  // namespace everybeam
#endif  // EVERYBEAM_CIRCULARSYMMETRIC_ATCACOEFFICIENTS_H_
