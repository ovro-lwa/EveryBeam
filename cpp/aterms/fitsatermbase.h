// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_FITS_ATERM_BASE_H
#define EVERYBEAM_ATERMS_FITS_ATERM_BASE_H

#include "atermbase.h"
#include "atermresampler.h"
#include "cache.h"

#include <aocommon/fits/fitsreader.h>
#include <aocommon/windowfunction.h>

namespace everybeam {
namespace aterms {

class FitsATermBase : public ATermBase {
 public:
  FitsATermBase(size_t n_antennas,
                const aocommon::CoordinateSystem& coordinate_system,
                size_t max_support);
  ~FitsATermBase() override;

  double AverageUpdateTime() const override;

  void SetTukeyWindow(double padding) { resampler_.SetTukeyWindow(padding); }

  void SetWindow(aocommon::WindowFunction::Type window) {
    resampler_.SetWindow(window);
  }

  void SetDownSample(bool downsample) { resampler_.SetDownSample(downsample); }

 protected:
  const ATermResampler& GetResampler() const { return resampler_; }

  void ReadAndResample(aocommon::FitsReader& reader, size_t file_index,
                       aocommon::UVector<float>& scratch,
                       aocommon::UVector<float>& output,
                       double stretch_factor) {
    resampler_.ReadAndResample(reader, file_index, scratch, output,
                               stretch_factor);
  }

  void InitializeFromFiles(std::vector<aocommon::FitsReader>& readers);

  bool FindFilePosition(std::complex<float>* buffer, double time,
                        double frequency, size_t& time_index,
                        bool& requires_recalculation);

  void StoreInCache(double frequency, const std::complex<float>* buffer);

  struct Timestep {
    Timestep(double t, size_t ri, size_t ii)
        : time(t), reader_index(ri), img_index(ii) {}
    double time;
    size_t reader_index;
    size_t img_index;
  };

  size_t Width() const { return coordinate_system_.width; }
  size_t Height() const { return coordinate_system_.height; }
  size_t NAntennas() const { return n_antennas_; }
  size_t NFrequencies() const { return n_frequencies_; }
  const Timestep& GetTimestep(size_t index) const { return timesteps_[index]; }
  const aocommon::CoordinateSystem GetCoordinateSystem() const {
    return coordinate_system_;
  }

 private:
  std::vector<Timestep> timesteps_;
  Cache cache_;
  size_t cur_timeindex_;
  double cur_frequency_;
  size_t n_frequencies_;
  const size_t n_antennas_;
  const aocommon::CoordinateSystem coordinate_system_;
  ATermResampler resampler_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
