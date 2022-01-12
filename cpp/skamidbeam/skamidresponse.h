// skamidresponse.h Base class for SKA-MID responses.
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_SKAMIDBEAM_SKAMIDRESPONSE_H_
#define EVERYBEAM_SKAMIDBEAM_SKAMIDRESPONSE_H_

#include <complex>

namespace everybeam {
namespace skamidbeam {
/**
 * @brief Base class for SKA-MID beam response models.
 *
 */
class SkaMidResponse {
 public:
  /**
   * @brief Construct a new Ska Mid Response object
   *
   * @param diameter Diameter of dish (m)
   * @param blockage Blockage of dish due to receiver (m)
   */
  SkaMidResponse(double diameter, double blockage)
      : diameter_(diameter), blockage_(blockage) {
    const double ratio = blockage_ / diameter_;
    blockage_factor_ = ratio * ratio;
  }

  virtual ~SkaMidResponse() = default;

  /**
   * @brief Render voltage pattern on grid
   *
   * @param aterm Buffer for storing the beam response. Should have size width x
   * height * 4.
   * @param width Width of grid.
   * @param height Height of grid.
   * @param pixel_scale_x Pixel scale in width direction (rad)
   * @param pixel_scale_y Pixel scale in height direction (rad)
   * @param phase_centre_ra Right ascension of phase centre (rad)
   * @param phase_centre_dec Declination of phase centre (rad)
   * @param pointing_ra Right ascension of telescope pointing (rad)
   * @param pointing_dec Declination of telescope pointing (rad)
   * @param phase_centre_dl Shift in l from pointing direction to grid centre
   * (rad)
   * @param phase_centre_dm Shift in m from pointing direction to grid centre
   * (rad)
   * @param frequency_hz Frequency (Hz)
   */
  virtual void Render(std::complex<float>* aterm, size_t width, size_t height,
                      double pixel_scale_x, double pixel_scale_y,
                      double phase_centre_ra, double phase_centre_dec,
                      double pointing_ra, double pointing_dec,
                      double phase_centre_dl, double phase_centre_dm,
                      double frequency_hz) const = 0;

  /**
   * @brief Render voltage pattern for single point.
   *
   * @param aterm Buffer for storing the beam response. Should accomodate 4
   * complex numbers (i.e. a Jones matrix).
   * @param phase_centre_ra Right ascension of phase centre (rad)
   * @param phase_centre_dec Declination of phase centre (rad)
   * @param pointing_ra Right ascension of telescope pointing (rad)
   * @param pointing_dec Declination of telescope pointing (rad)
   * @param frequency_hz Frequency (Hz)
   */
  virtual void Render(std::complex<float>* aterm, double phase_centre_ra,
                      double phase_centre_dec, double pointing_ra,
                      double pointing_dec, double frequency_hz) const = 0;

 protected:
  double diameter_;
  double blockage_;
  double blockage_factor_;
};
}  // namespace skamidbeam
}  // namespace everybeam

#endif