// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_DLDM_ATERM_H
#define EVERYBEAM_DLDM_ATERM_H

#include "fitsatermbase.h"

#include <aocommon/fits/fitsreader.h>

namespace everybeam {
namespace aterms {

class DLDMATerm final : public FitsATermBase {
 public:
  DLDMATerm(size_t n_antenna,
            const aocommon::CoordinateSystem& coordinate_system,
            size_t max_support);

  void Open(const std::vector<std::string>& filenames);

  bool Calculate(std::complex<float>* buffer, double time, double frequency,
                 size_t fieldId, const double* uvwInM) override;

  void SetUpdateInterval(double updateInterval) {
    update_interval_ = updateInterval;
  }

  double AverageUpdateTime() const final override {
    return std::min(FitsATermBase::AverageUpdateTime(), update_interval_);
  }

 private:
  void ReadImages(std::complex<float>* buffer, size_t timeIndex,
                  double frequency, const double* uvwInM);
  void EvaluateDLDM(std::complex<float>* dest, const float* dl, const float* dm,
                    const double* uvwInM);

  std::vector<aocommon::FitsReader> readers_;
  aocommon::UVector<float> scratch_;
  aocommon::UVector<float> dl_image_;
  aocommon::UVector<float> dm_image_;
  std::vector<std::array<double, 2>> uvws_;
  double update_interval_;
  double previous_time_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
