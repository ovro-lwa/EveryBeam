// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_FITS_ATERM_H
#define EVERYBEAM_ATERMS_FITS_ATERM_H

#include "fitsatermbase.h"

#include <aocommon/fits/fitsreader.h>
#include <aocommon/uvector.h>

#include <complex>
#include <map>
#include <vector>

namespace everybeam {
namespace aterms {

/**
 * Class that reads in FITS images and resamples them onto aterm grids.
 * The fits file is supposed to have a TIME, FREQ and ANTENNA axis.
 */
class FitsATerm final : public FitsATermBase {
 public:
  FitsATerm(size_t n_antennas,
            const aocommon::CoordinateSystem& coordinate_system,
            size_t max_support);
  ~FitsATerm() override;

  void OpenTECFiles(const std::vector<std::string>& filenames);
  void OpenDiagGainFiles(const std::vector<std::string>& filenames);

  bool Calculate(std::complex<float>* buffer, double time, double frequency,
                 size_t field_id, const double* uvw_in_m) override;

 private:
  enum class Mode { kTEC, kDiagonal } mode_;

  void ReadImages(std::complex<float>* buffer, size_t time_index,
                  double frequency);

  // TODO: Either implement and use or remove this function.
  // void Resample(const aocommon::FitsReader& reader, double* dest,
  //              const double* source);

  void EvaluateTEC(std::complex<float>* dest, const float* source,
                   double frequency);

  void CopyToRealPolarization(std::complex<float>* dest, const float* source,
                              size_t pol_index);
  void CopyToImaginaryPolarization(std::complex<float>* dest,
                                   const float* source, size_t pol_index);
  void SetPolarization(std::complex<float>* dest, size_t pol_index,
                       std::complex<float> value);

  aocommon::UVector<float> scratch_a_, scratch_b_;
  std::vector<aocommon::FitsReader> readers_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
