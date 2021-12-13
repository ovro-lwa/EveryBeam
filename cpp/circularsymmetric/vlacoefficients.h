// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_CIRCULARSYMMETRIC_VLACOEFFICIENTS_H_
#define EVERYBEAM_CIRCULARSYMMETRIC_VLACOEFFICIENTS_H_

#include <array>
#include <map>
#include <string>

#include "coefficients.h"

namespace everybeam {
namespace circularsymmetric {
class VLACoefficients final : public Coefficients {
 public:
  VLACoefficients(const std::string& band_name) : band_name_(band_name) {}

  aocommon::UVector<double> GetFrequencies(double frequency) const override {
    return aocommon::UVector<double>{frequency};
  }
  aocommon::UVector<double> GetCoefficients(double frequency) const override {
    std::array<double, 5> coefficients = GetCoefficients(band_name_, frequency);
    return aocommon::UVector<double>(coefficients.begin(), coefficients.end());
  }
  double MaxRadiusInArcMin() const override { return 53.0; }
  double ReferenceFrequency() const override { return 1e9; }
  bool AreInverted() const override { return false; }

 private:
  const std::string band_name_;

  static std::array<double, 5> GetCoefficients(const std::string& band_name,
                                               double freq);
  static std::map<int, std::array<double, 5>> GetCoefficients();
  static std::map<char, double> GetFeedConf();
  static char DetermineFeed(double freq, double freq_center = 0.0);
  static void LimitFreqForBand(char band, double& freq);
};
}  // namespace circularsymmetric
}  // namespace everybeam
#endif  // EVERYBEAM_CIRCULARSYMMETRIC_VLACOEFFICIENTS_H_
