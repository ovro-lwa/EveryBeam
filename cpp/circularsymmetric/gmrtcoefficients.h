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
 * - https://github.com/ruta-k/uGMRTprimarybeam
 * - PBCOR in AIPS and code/synthesis/TransformMachines/PBMath.cc in CASA.
 */
class GMRTCoefficients final : public Coefficients {
 public:
  aocommon::UVector<double> GetFrequencies(
      [[maybe_unused]] double frequency) const override {
    return {375e6, 700e6, 1250e6};
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
  // These come from
  // "Upgraded GMRT: Preliminary Primary-beam shape parameters",
  // http://www.ncra.tifr.res.in/ncra/gmrt/gmrt-users/observing-help/ugmrt-primary-beam-shape
  // coef x nfreq
  static constexpr std::array<double, 15> coefficients_{
      1.0, -2.939e-03, 33.312e-07, -16.659e-10, 3.066e-13,
      1.0, -3.190e-03, 38.642e-07, -20.471e-10, 3.964e-13,
      1.0, -2.608e-03, 27.357e-07, -13.091e-10, 2.368e-13,
  };
};
}  // namespace circularsymmetric
}  // namespace everybeam
#endif  // EVERYBEAM_CIRCULARSYMMETRIC_ATCACOEFFICIENTS_H_
