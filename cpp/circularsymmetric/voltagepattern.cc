// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "voltagepattern.h"

#include <aocommon/imagecoordinates.h>

#include <cmath>

using aocommon::ImageCoordinates;
using aocommon::UVector;
using everybeam::circularsymmetric::VoltagePattern;

void VoltagePattern::EvaluatePolynomial(const UVector<double>& coefficients,
                                        bool invert) {
  // This comes from casa's: void PBMath1DIPoly::fillPBArray(), wideband case
  size_t nsamples = 10000;
  size_t nfreq = frequencies_.size();
  size_t ncoef = coefficients.size() / nfreq;
  values_.resize(nsamples * nfreq);
  inverse_increment_radius_ = double(nsamples - 1) / maximum_radius_arc_min_;
  double* output = values_.data();
  for (size_t n = 0; n != nfreq; n++) {
    const double* freqcoefficients = &coefficients[n * ncoef];
    for (size_t i = 0; i < nsamples; i++) {
      double taper = 0.0;
      double x2 = double(i) / inverse_increment_radius_;
      x2 = x2 * x2;
      double y = 1.0;

      for (size_t j = 0; j < ncoef; j++) {
        taper += y * freqcoefficients[j];
        y *= x2;
      }
      if (taper >= 0.0) {
        if (invert && taper != 0.0) {
          taper = 1.0 / std::sqrt(taper);
        } else {
          taper = std::sqrt(taper);
        }
      } else {
        taper = 0.0;
      }
      *output = taper;
      ++output;
    }
  }
}

UVector<double> VoltagePattern::InterpolateValues(double freq) const {
  UVector<double> result;
  size_t ifit = 0;
  size_t nfreq = frequencies_.size();
  for (ifit = 0; ifit != nfreq; ifit++) {
    if (freq <= frequencies_[ifit]) break;
  }
  if (ifit == 0) {
    result.assign(values_.begin(), values_.begin() + NSamples());
  } else if (ifit == nfreq) {
    result.assign(values_.begin() + (nfreq - 1) * NSamples(), values_.end());
  } else {
    size_t n = NSamples();
    double l = (freq - frequencies_[ifit - 1]) /
               (frequencies_[ifit] - frequencies_[ifit - 1]);
    const double* vpA = FreqIndexValues(ifit - 1);
    const double* vpB = FreqIndexValues(ifit);
    result.resize(n);
    for (size_t i = 0; i != n; ++i) {
      result[i] = vpA[i] * (1.0 - l) + vpB[i] * l;
    }
  }
  return result;
}

const double* VoltagePattern::InterpolateValues(
    double frequency_hz, UVector<double>& interpolated_values) const {
  if (frequencies_.size() > 1) {
    interpolated_values = InterpolateValues(frequency_hz);
    return interpolated_values.data();
  } else {
    return FreqIndexValues(0);
  }
}

double VoltagePattern::LmMaxSquared(double frequency_hz) const {
  double factor =
      (180.0 / M_PI) * 60.0 * frequency_hz * 1.0e-9;  // arcminutes * GHz
  double rmax = maximum_radius_arc_min_ / factor;
  return rmax * rmax;
}

void VoltagePattern::Render(std::complex<float>* aterm, size_t width,
                            size_t height, double pixel_scale_x,
                            double pixel_scale_y, double phase_centre_ra,
                            double phase_centre_dec, double pointing_ra,
                            double pointing_dec, double phase_centre_dl,
                            double phase_centre_dm, double frequency_hz) const {
  double lmMaxSq = LmMaxSquared(frequency_hz);

  UVector<double> interpolated_values;
  const double* vp = InterpolateValues(frequency_hz, interpolated_values);

  double factor =
      (180.0 / M_PI) * 60.0 * frequency_hz * 1.0e-9;  // arcminutes * GHz
  double l0, m0;
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
      l += phase_centre_dl;
      m += m0;
      ImageCoordinates::LMToRaDec(l, m, phase_centre_ra, phase_centre_dec, ra,
                                  dec);
      ImageCoordinates::RaDecToLM(ra, dec, pointing_ra, pointing_dec, l, m);
      l -= l0;
      m -= m0;
      double r2 = l * l + m * m;
      double out;
      if (r2 > lmMaxSq) {
        out = 1e-4;
      } else {
        double r = std::sqrt(r2) * factor;
        int indx = int(r * inverse_increment_radius_);
        out = vp[indx] * (1.0 - 1e-4) + 1e-4;
      }

      std::complex<float>* ptr = row + ix * 4;
      ptr[0] = out;
      ptr[1] = 0.0;
      ptr[2] = 0.0;
      ptr[3] = out;
    }
  }
}

void VoltagePattern::Render(std::complex<float>* aterm, double phase_centre_ra,
                            double phase_centre_dec, double pointing_ra,
                            double pointing_dec, double frequency_hz) const {
  double lmMaxSq = LmMaxSquared(frequency_hz);

  UVector<double> interpolated_values;
  const double* vp = InterpolateValues(frequency_hz, interpolated_values);

  double factor =
      (180.0 / M_PI) * 60.0 * frequency_hz * 1.0e-9;  // arcminutes * GHz

  // TODO: probably not all conversions needed?
  double l0, m0;
  ImageCoordinates::RaDecToLM(pointing_ra, pointing_dec, phase_centre_ra,
                              phase_centre_dec, l0, m0);
  double l = 0;
  double m = m0;
  double ra, dec;
  ImageCoordinates::LMToRaDec(l, m, phase_centre_ra, phase_centre_dec, ra, dec);
  ImageCoordinates::RaDecToLM(ra, dec, pointing_ra, pointing_dec, l, m);
  l -= l0;
  m -= m0;
  double r2 = l * l + m * m;
  double out;
  if (r2 > lmMaxSq) {
    out = 1e-4;
  } else {
    double r = std::sqrt(r2) * factor;
    int indx = int(r * inverse_increment_radius_);
    out = vp[indx] * (1.0 - 1e-4) + 1e-4;
  }

  aterm[0] = out;
  aterm[1] = 0.0;
  aterm[2] = 0.0;
  aterm[3] = out;
}
