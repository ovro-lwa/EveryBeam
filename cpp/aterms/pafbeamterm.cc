// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "pafbeamterm.h"

#include <boost/algorithm/string.hpp>

#include <aocommon/fits/fitsreader.h>

namespace everybeam {
namespace aterms {

PAFBeamTerm::PAFBeamTerm(const aocommon::CoordinateSystem& coordinate_system,
                         size_t max_support)
    : coordinate_system_(coordinate_system),
      resampler_(coordinate_system, max_support),
      freq0_(0.0),
      dfreq_(0.0),
      beam_ra_(0.0),
      beam_dec_(0.0),
      update_interval_(3600),
      previous_time_(std::numeric_limits<double>::lowest()),
      correct_for_frequency_offset_(true),
      ref_frequency_(0.0) {}

bool PAFBeamTerm::Calculate(std::complex<float>* buffer, double time,
                            double frequency, size_t /*field_id*/,
                            const double* /*uvw_in_m*/) {
  bool outdated = std::fabs(time - previous_time_) > update_interval_;
  if (!outdated) return false;
  previous_time_ = time;

  const size_t width = coordinate_system_.width;
  const size_t height = coordinate_system_.height;
  const size_t freq_index = std::min<size_t>(
      n_frequencies_ - 1, std::max(0.0, round((frequency - freq0_) / dfreq_)));

  const double ref_frequency =
      (ref_frequency_ != 0.0) ? ref_frequency_ : (freq0_ + freq_index * dfreq_);
  const double frequency_factor =
      correct_for_frequency_offset_ ? frequency / ref_frequency : 1.0;

  aocommon::UVector<float> scratch(resampler_.ScratchASize());
  aocommon::UVector<float> output(resampler_.ScratchBSize(readers_[0]));
  for (size_t antenna = 0; antenna != n_antennas_; ++antenna) {
    resampler_.OverrideFitsPhaseCentre(beam_ra_, beam_dec_);
    resampler_.ReadAndResample(readers_[antenna], freq_index, scratch, output,
                               frequency_factor);
    for (size_t i = 0; i != width * height; ++i) {
      buffer[0] = output[i];
      buffer[1] = 0.0;
      buffer[2] = 0.0;
      buffer[3] = output[i];
      buffer += 4;
    }
  }

  return true;
}

void PAFBeamTerm::Open(const std::string& filename_template,
                       const std::vector<std::string>& antenna_map,
                       const std::string& beam_name, double beam_ra,
                       double beam_dec) {
  n_antennas_ = antenna_map.size();
  readers_.clear();
  beam_ra_ = beam_ra;
  beam_dec_ = beam_dec;
  for (size_t antenna = 0; antenna != n_antennas_; ++antenna) {
    std::string filename = boost::algorithm::replace_all_copy(
        filename_template, "$ANT", antenna_map[antenna]);
    boost::algorithm::replace_all(filename, "$BEAM", beam_name);
    readers_.emplace_back(filename, false, true);
    if (antenna == 0) {
      freq0_ = readers_.back().FrequencyDimensionStart();
      dfreq_ = readers_.back().FrequencyDimensionIncr();
      n_frequencies_ = readers_.back().NFrequencies();
    } else {
      if (freq0_ != readers_.back().FrequencyDimensionStart() ||
          dfreq_ != readers_.back().FrequencyDimensionIncr() ||
          n_frequencies_ != readers_.back().NFrequencies()) {
        throw std::runtime_error(
            "All antenna fits files should have the same frequency axis");
      }
    }
  }
}

}  // namespace aterms
}  // namespace everybeam
