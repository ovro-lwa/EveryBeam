// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_SKAMIDBEAM_SKAMIDANALYTICALRESPONSE_H_
#define EVERYBEAM_SKAMIDBEAM_SKAMIDANALYTICALRESPONSE_H_

#include "skamidresponse.h"

#include <complex>

namespace everybeam {
namespace skamidbeam {
/**
 * @brief Analtical and circularsymmetric beam response model
 * for a SKA-MID dish.
 *
 */
class SkaMidAnalyticalResponse final : public SkaMidResponse {
 public:
  [[gnu::visibility("default")]] SkaMidAnalyticalResponse(double diameter,
                                                          double blockage);

  void Render(std::complex<float>* aterm, size_t width, size_t height,
              double pixel_scale_x, double pixel_scale_y,
              double phase_centre_ra, double phase_centre_dec,
              double pointing_ra, double pointing_dec, double phase_centre_dl,
              double phase_centre_dm, double frequency_hz) const override;

  void Render(std::complex<float>* aterm, double phase_centre_ra,
              double phase_centre_dec, double pointing_ra, double pointing_dec,
              double frequency_hz) const override;

 private:
  /**
   * @brief Voltage pattern is based on evaluating a Bessel function of the
   * first kind.
   *
   * @param r Radius (m)
   * @return float
   */
  float EvaluateBessel(double r) const;
};
}  // namespace skamidbeam
}  // namespace everybeam
#endif