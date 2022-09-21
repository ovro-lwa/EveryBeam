// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_ATERM_RESAMPLER_H
#define EVERYBEAM_ATERMS_ATERM_RESAMPLER_H

#include "atermbase.h"

#include <aocommon/coordinatesystem.h>
#include <aocommon/uvector.h>
#include <aocommon/windowfunction.h>

namespace aocommon {
class FitsReader;
}

namespace everybeam {

namespace common {
class FFTResampler;
}

namespace aterms {

class ATermResampler {
 public:
  ATermResampler(const aocommon::CoordinateSystem& coordinate_system,
                 size_t max_support);
  ~ATermResampler();

  /**
   * @param scratch vector of size at least ScratchASize()
   * @param output vector of size at least ScratchBSize()
   */
  void ReadAndResample(aocommon::FitsReader& reader, size_t file_index,
                       aocommon::UVector<float>& scratch,
                       aocommon::UVector<float>& output, double stretch_factor);

  void SetTukeyWindow(double padding) {
    window_ = aocommon::WindowFunction::Tukey;
    padding_ = padding;
  }

  void SetWindow(aocommon::WindowFunction::Type window) { window_ = window; }

  void SetDownSample(bool downsample) { downsample_ = downsample; }

  size_t AllocatedWidth() const { return allocated_width_; }
  size_t AllocatedHeight() const { return allocated_height_; }

  size_t ScratchASize() const { return allocated_width_ * allocated_height_; }

  size_t ScratchBSize(const aocommon::FitsReader& reader) const;

  void OverrideFitsPhaseCentre(double ra, double dec) {
    override_fits_phase_centre_ = true;
    override_ra_ = ra;
    override_dec_ = dec;
  }

 private:
  void regrid(const aocommon::FitsReader& reader, float* dest,
              const float* source, double stretch_factor);

  const aocommon::CoordinateSystem coordinate_system_;
  size_t allocated_width_;
  size_t allocated_height_;
  std::unique_ptr<common::FFTResampler> resampler_;
  bool downsample_;
  aocommon::WindowFunction::Type window_;
  double padding_;
  bool override_fits_phase_centre_;
  double override_ra_;
  double override_dec_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
