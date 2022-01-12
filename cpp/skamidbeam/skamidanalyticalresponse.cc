// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "skamidanalyticalresponse.h"
#include "../common/constants.h"

#include <aocommon/imagecoordinates.h>

#include <cmath>

#include <boost/math/special_functions/bessel.hpp>

using aocommon::ImageCoordinates;

namespace everybeam {
namespace skamidbeam {
SkaMidAnalyticalResponse::SkaMidAnalyticalResponse(double diameter,
                                                   double blockage)
    : SkaMidResponse(diameter, blockage) {}

void SkaMidAnalyticalResponse::Render(
    std::complex<float>* aterm, size_t width, size_t height,
    double pixel_scale_x, double pixel_scale_y, double phase_centre_ra,
    double phase_centre_dec, double pointing_ra, double pointing_dec,
    double phase_centre_dl, double phase_centre_dm, double frequency_hz) const {
  const double wavelength = common::c / frequency_hz;

  double l0;
  double m0;
  ImageCoordinates::RaDecToLM(pointing_ra, pointing_dec, phase_centre_ra,
                              phase_centre_dec, l0, m0);
  l0 += phase_centre_dl;
  m0 += phase_centre_dm;

  for (size_t iy = 0; iy != height; ++iy) {
    std::complex<float>* row = aterm + iy * width * 4;
    for (size_t ix = 0; ix != width; ++ix) {
      double l, m, ra, dec;
      ImageCoordinates::XYToLM(ix, iy, pixel_scale_x, pixel_scale_y, width,
                               height, l, m);
      l += l0;
      m += m0;
      ImageCoordinates::LMToRaDec(l, m, phase_centre_ra, phase_centre_dec, ra,
                                  dec);
      ImageCoordinates::RaDecToLM(ra, dec, pointing_ra, pointing_dec, l, m);
      l -= l0;
      m -= m0;

      // TODO: no check on lm_max_sq needed?
      const double r = std::sqrt(l * l + m * m);
      const float vp_reflector =
          EvaluateBessel(r * M_PI * diameter_ / wavelength);
      const float vp_blockage =
          EvaluateBessel(r * M_PI * blockage_ / wavelength);
      // NOTE: in rascil create_pb() the voltage pattern is multiplied with its
      // complex conjugate. In WSClean - EveryBeam, this multiplication is done
      // by the callee.
      const std::complex<float> vp_combined =
          vp_reflector - blockage_factor_ * vp_blockage;

      std::complex<float>* ptr = row + ix * 4;
      ptr[0] = vp_combined;
      ptr[1] = 0.0;
      ptr[2] = 0.0;
      ptr[3] = vp_combined;
    }
  }
}

void SkaMidAnalyticalResponse::Render(std::complex<float>* aterm,
                                      double phase_centre_ra,
                                      double phase_centre_dec,
                                      double pointing_ra, double pointing_dec,
                                      double frequency_hz) const {
  const double wavelength = common::c / frequency_hz;

  double l0;
  double m0;
  ImageCoordinates::RaDecToLM(pointing_ra, pointing_dec, phase_centre_ra,
                              phase_centre_dec, l0, m0);
  double l = l0;
  double m = m0;
  double ra, dec;
  ImageCoordinates::LMToRaDec(l, m, phase_centre_ra, phase_centre_dec, ra, dec);
  ImageCoordinates::RaDecToLM(ra, dec, pointing_ra, pointing_dec, l, m);

  l -= l0;
  m -= m0;
  // TODO: no check on lm_max_sq needed?
  const double r = std::sqrt(l * l + m * m);
  const float vp_reflector = EvaluateBessel(r * M_PI * diameter_ / wavelength);
  const float vp_blockage = EvaluateBessel(r * M_PI * blockage_ / wavelength);
  // NOTE: in rascil create_pb() the voltage pattern is multiplied with its
  // complex conjugate. In WSClean - EveryBeam, this multiplication is done by
  // the callee.
  const std::complex<float> vp_combined =
      vp_reflector - blockage_factor_ * vp_blockage;

  // TODO: different polarisations?
  aterm[0] = vp_combined;
  aterm[1] = 0.0;
  aterm[2] = 0.0;
  aterm[3] = vp_combined;
}

// std::complex<double> ?
float SkaMidAnalyticalResponse::EvaluateBessel(double r) const {
  const float order = 1.0;
  const double r_small = 1e-9;
  const double radius = (r > 0.0) ? r : r_small;
  return 2.0 * boost::math::cyl_bessel_j(order, radius) / radius;
}
}  // namespace skamidbeam
}  // namespace everybeam