// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_PAF_BEAM_TERM_H
#define EVERYBEAM_ATERMS_PAF_BEAM_TERM_H

#include "atermbase.h"
#include "atermresampler.h"

#include <aocommon/fits/fitsreader.h>

#include <string>
#include <vector>

namespace everybeam {
namespace aterms {

/**
 * Class for reading fits files as they are used in "Phased-array feed"
 * such as Apertif.
 */
class PAFBeamTerm final : public ATermBase {
 public:
  PAFBeamTerm(const aocommon::CoordinateSystem& coordinate_system,
              size_t max_support);

  bool Calculate(std::complex<float>* buffer, double time, double frequency,
                 size_t fieldId, const double* uvw_in_m) override;

  void SetUpdateInterval(double update_interval) {
    update_interval_ = update_interval;
  }

  double AverageUpdateTime() const final override { return update_interval_; }

  void Open(const std::string& filename_template,
            const std::vector<std::string>& antenna_map,
            const std::string& beam_name, double beam_ra, double beam_dec);

  void SetTukeyWindow(double padding) { resampler_.SetTukeyWindow(padding); }

  void SetWindow(aocommon::WindowFunction::Type window) {
    resampler_.SetWindow(window);
  }

  void SetDownSample(bool downsample) { resampler_.SetDownSample(downsample); }

  void SetCorrectForFrequencyOffset(bool correct) {
    correct_for_frequency_offset_ = correct;
  }
  void SetReferenceFrequency(double ref_frequency) {
    ref_frequency_ = ref_frequency;
  }

 private:
  std::vector<aocommon::FitsReader> readers_;
  const aocommon::CoordinateSystem coordinate_system_;
  ATermResampler resampler_;
  size_t n_antennas_;
  size_t n_frequencies_;
  double freq0_;
  double dfreq_;
  double beam_ra_;
  double beam_dec_;
  double update_interval_;
  double previous_time_;
  bool correct_for_frequency_offset_;
  double ref_frequency_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
