// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "config.h"
#include "../aterms/h5parmaterm.h"
#include "../coords/coordutils.h"

#include <aocommon/uvector.h>
#include <aocommon/imagecoordinates.h>
#include <vector>
#include <complex>
#include <math.h>

using aocommon::CoordinateSystem;
using everybeam::aterms::H5ParmATerm;
using everybeam::aterms::LagrangePolynomial;
using schaapcommon::h5parm::H5Parm;

std::string h5parm_mock = MOCK_H5PARM_PATH;

BOOST_AUTO_TEST_SUITE(th5parmaterm)
// Explicit, but inefficient evaluation of polynomials, for reference only
float ExplicitPolynomialEvaluation(float x, float y, size_t order,
                                   const std::vector<float>& coeffs) {
  float sol = 0;
  size_t idx = 0;
  for (size_t n = 0; n < order + 1; ++n) {
    for (size_t k = 0; k < n + 1; ++k) {
      sol += coeffs[idx] * pow(x, n - k) * pow(y, k);
      idx += 1;
    }
  }
  return sol;
}

// Convenience method to convert pixel coordinates to lm coordinates
std::vector<std::pair<float, float>> ConvertXYToLM(
    CoordinateSystem coord_system) {
  std::vector<std::pair<float, float>> result;
  for (size_t y = 0; y < coord_system.height; ++y) {
    for (size_t x = 0; x < coord_system.width; ++x) {
      double l, m;
      aocommon::ImageCoordinates::XYToLM(x, y, coord_system.dl, coord_system.dm,
                                         coord_system.width,
                                         coord_system.height, l, m);
      result.push_back(std::make_pair(l, m));
    }
  }
  return result;
}

/**
 * @brief Compute coefficients that should be stored in the MOCK_H5PARM.h5 file.
 * Compute coefficients that should be stored in the MOCK_H5PARM.h5 file.
 * flattened, where TOTAL_NR_COEFFS = nstations * ntimes * ncoeffs
 *
 * The (python) script for generating the MOCK_H5PARM.h5 file can be found in
 * the scripts/misc/ directory (scripts/misc/make_mock_h5parm.py).
 *
 * @param porder Polynomial order
 * @param sidx Station index
 * @param tidx Time index
 * @return std::vector<float> Vector of coefficients
 */
std::vector<float> ComputeTestCoeffs(size_t porder, size_t sidx, hsize_t tidx) {
  size_t nr_coeffs = LagrangePolynomial::ComputeNrCoeffs(porder);
  LagrangePolynomial polynomial(nr_coeffs);
  std::vector<float> coeffs(nr_coeffs);
  for (size_t i = 0; i < coeffs.size(); ++i) {
    coeffs[i] = (nr_coeffs * sidx + (i + 1)) * (tidx + 1);
  }
  return coeffs;
}

// Zero order polynomial
BOOST_AUTO_TEST_CASE(test_lagrange_0) {
  LagrangePolynomial poly(1);
  BOOST_CHECK_EQUAL(poly.GetOrder(), 0);

  // Make coefficients for binomial:
  // f(x,y) = 10
  std::vector<float> coeffs = {10};

  float x = 20, y = 6;
  float result = poly.Evaluate(x, y, coeffs);
  float ref = ExplicitPolynomialEvaluation(x, y, poly.GetOrder(), coeffs);
  BOOST_CHECK_EQUAL(result, ref);
}

// First order polynomial
BOOST_AUTO_TEST_CASE(test_lagrange_1) {
  LagrangePolynomial poly(3);
  BOOST_CHECK_EQUAL(poly.GetOrder(), 1);

  // Make coefficients for binomial:
  // f(x,y) = 10 - 4x + 3y
  std::vector<float> coeffs = {10, -4, 3};

  float x = 20, y = 6;
  float result = poly.Evaluate(x, y, coeffs);
  float ref = ExplicitPolynomialEvaluation(x, y, poly.GetOrder(), coeffs);
  BOOST_CHECK_EQUAL(result, ref);
}

// Second order polynomial
BOOST_AUTO_TEST_CASE(test_lagrange_2) {
  // Second order binomial has 6 terms
  LagrangePolynomial poly(6);
  BOOST_CHECK_EQUAL(poly.GetOrder(), 2);

  // Make coefficients for binomial:
  // f(x,y) = 4 + 3x - 2y + 5x^2 - 2xy + 4y^2
  std::vector<float> coeffs = {4, 3, -2, 5, -2, 4};

  float x = 20, y = 6;
  float result = poly.Evaluate(x, y, coeffs);
  float ref = ExplicitPolynomialEvaluation(x, y, poly.GetOrder(), coeffs);
  BOOST_CHECK_EQUAL(result, ref);
}

// Third order polynomial
BOOST_AUTO_TEST_CASE(test_lagrange_3) {
  // Third order binomial has 10 terms
  LagrangePolynomial poly(10);
  BOOST_CHECK_EQUAL(poly.GetOrder(), 3);

  // Make coefficients for binomial:
  // f(x,y) = 10 - x + 2y + 3x^2 + 4xy - 5y^2 + 6x^3 + 0x^2y + 7xy^2 + 8y^3
  std::vector<float> coeffs = {10, -1, 2, 3, 4, -5, 6, 2, 7, 8};

  float x = 20, y = 6;
  float result = poly.Evaluate(x, y, coeffs);
  float ref = ExplicitPolynomialEvaluation(x, y, poly.GetOrder(), coeffs);
  BOOST_CHECK_EQUAL(result, ref);
}

BOOST_AUTO_TEST_CASE(read_h5parmfile) {
  std::vector<std::string> h5parm_files = {h5parm_mock};
  // Names in MOCK_H5PARM.h5 should match Antenna0, Antenna1
  std::vector<std::string> station_names = {"Antenna0", "Antenna1"};

  double frequency = 57812500.;
  double ra(-1.44194878), dec(0.85078091);

  // Properties of grid
  std::size_t width(4), height(4);
  double dl(0.5 * M_PI / 180.), dm(0.5 * M_PI / 180.), shift_l(0.), shift_m(0.);

  CoordinateSystem coord_system = {width, height, ra,      dec,
                                   dl,    dm,     shift_l, shift_m};

  H5ParmATerm h5parmaterm(station_names, coord_system);
  h5parmaterm.Open(h5parm_files);

  // Needed for reference solution
  H5Parm h5parm_tmp = H5Parm(h5parm_files[0]);
  LagrangePolynomial ampl_polynomial(6);
  LagrangePolynomial phase_polynomial(3);
  std::vector<std::pair<float, float>> image_coords =
      ConvertXYToLM(coord_system);

  for (float time = 0; time <= 10; ++time) {
    // Compute solution with H5ParmATerm
    aocommon::UVector<std::complex<float>> h5parm_buffer(
        width * height * station_names.size() * 4);
    h5parmaterm.Calculate(h5parm_buffer.data(), time, frequency, 0, nullptr);

    // Compute reference solution
    hsize_t tindex_ampl =
        h5parm_tmp.GetSolTab("amplitude_coefficients").GetTimeIndex(time);
    hsize_t tindex_phase =
        h5parm_tmp.GetSolTab("phase_coefficients").GetTimeIndex(time);

    for (size_t i = 0; i < station_names.size(); ++i) {
      std::vector<float> ampl_coeffs = ComputeTestCoeffs(2, i, tindex_ampl);
      std::vector<float> phase_coeffs = ComputeTestCoeffs(1, i, tindex_phase);
      for (size_t j = 0; j < width * height; ++j) {
        const float l = image_coords[j].first;
        const float m = image_coords[j].second;

        float ampl_ref = ampl_polynomial.Evaluate(l, m, ampl_coeffs);
        float phase_ref = phase_polynomial.Evaluate(l, m, phase_coeffs);
        std::complex<float> ref =
            ampl_ref * exp(std::complex<float>(0, phase_ref));

        size_t offset = (i * width * height + j) * 4;
        std::complex<float> result = h5parm_buffer[offset];

        BOOST_CHECK_CLOSE(std::abs(result), std::abs(ref), 1e-4);
        BOOST_CHECK_CLOSE(std::arg(result), std::arg(ref), 1e-4);
      }
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()
