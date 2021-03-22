// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_ATERM_BEAM_H
#define EVERYBEAM_ATERMS_ATERM_BEAM_H

#include "atermbase.h"

namespace everybeam {
namespace aterms {
class ATermBeam : public ATermBase {
 public:
  ATermBeam()
      : update_interval_(0),
        last_aterm_update_(0),
        last_frequency_(0.0),
        last_field_id_(0) {}

  bool Calculate(std::complex<float>* buffer, double time, double frequency,
                 size_t field_id, const double*) final override {
    if (time - last_aterm_update_ > update_interval_ ||
        field_id != last_field_id_ || frequency != last_frequency_) {
      last_aterm_update_ = time;
      last_field_id_ = field_id;
      last_frequency_ = frequency;
      return CalculateBeam(buffer, time + update_interval_ * 0.5, frequency,
                           field_id);
    } else {
      return false;
    }
  }

  void SetUpdateInterval(double update_interval) {
    update_interval_ = update_interval;
    last_aterm_update_ = -update_interval_ - 1;
    last_field_id_ = 0;
  }

  double AverageUpdateTime() const final override { return update_interval_; }

 protected:
  virtual bool CalculateBeam(std::complex<float>* buffer, double time,
                             double frequency, size_t field_id) = 0;

 private:
  double update_interval_;
  double last_aterm_update_;
  double last_frequency_;
  size_t last_field_id_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
