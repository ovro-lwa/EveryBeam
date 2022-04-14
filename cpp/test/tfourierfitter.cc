// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "config.h"
#include "../aterms/fourierfitter.h"
#include <xtensor/xio.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor-fftw/basic.hpp>

#include <vector>
#include <complex>
#include <math.h>

using everybeam::aterms::FourierFitter;

#define TOLERANCE 1e-5

BOOST_AUTO_TEST_SUITE(tfourierfitter)

BOOST_AUTO_TEST_CASE(calculate_int_dirs) {
  const int kSubgridSize = 8;
  const int kSupport = 3;

  // Set of directiontions at exact integer pixels
  std::vector<std::pair<int, int>> directions_int = {{0, 0}, {0, 3}, {0, 6}};

  std::vector<std::pair<float, float>> directions(directions_int.begin(),
                                                  directions_int.end());

  std::vector<size_t> shape = {kSubgridSize, kSubgridSize};
  xt::xarray<std::complex<float>, xt::layout_type::row_major> screen(shape);

  FourierFitter fourier_fitter(kSubgridSize, kSupport, directions);

  const std::vector<std::complex<float>> solutions = {1.0, 2.0, 3.0};

  fourier_fitter.Evaluate(solutions, screen.data());

  for (size_t k = 0; k < directions.size(); ++k) {
    auto [i, j] = directions_int[k];
    BOOST_CHECK(abs(screen(i, j).real() - solutions[k].real()) < TOLERANCE);
    BOOST_CHECK(abs(screen(i, j).imag() - solutions[k].imag()) < TOLERANCE);
  }
}

BOOST_AUTO_TEST_CASE(calculate_float_dirs) {
  const int kSubgridSize = 8;
  const int kSupport = 3;
  const int kOversampling = 10;

  // Set of directiontions at as integer pixels positions in oversampled screen
  std::vector<std::pair<int, int>> directions_oversampled = {
      {5, 8}, {7, 32}, {3, 61}};

  // Transform integer positions in oversampled screen to floating point
  // positions in the (not oversampled) screen
  std::vector<std::pair<float, float>> directions;
  std::transform(directions_oversampled.begin(), directions_oversampled.end(),
                 std::back_inserter(directions),
                 [](std::pair<int, int> d) -> std::pair<float, float> {
                   return std::make_pair(float(d.first) / kOversampling,
                                         float(d.second) / kOversampling);
                 });

  std::vector<size_t> shape = {kSubgridSize, kSubgridSize};
  xt::xarray<std::complex<float>, xt::layout_type::row_major> screen(shape);

  FourierFitter fourier_fitter(kSubgridSize, kSupport, directions);

  const std::vector<std::complex<float>> solutions = {1.0, 2.0, 3.0};

  fourier_fitter.Evaluate(solutions, screen.data());

  // Oversample screen by kOversampling factor

  auto screen_ft = xt::fftw::fft2(screen);
  xt::xarray<std::complex<float>> screen_ft_oversampled =
      xt::zeros<std::complex<float>>(
          {kSubgridSize * kOversampling, kSubgridSize * kOversampling});
  for (int i = -kSubgridSize / 2; i < kSubgridSize / 2; ++i) {
    for (int j = -kSubgridSize / 2; j < kSubgridSize / 2; ++j) {
      screen_ft_oversampled.periodic(i, j) =
          screen_ft.periodic(i, j) *
          std::complex<float>(kOversampling * kOversampling);
    }
  }
  auto screen_oversampled = xt::fftw::ifft2(screen_ft_oversampled);

  // Compare values in oversampled screen to the requested value at given
  // directions

  for (size_t k = 0; k < directions.size(); ++k) {
    auto [i, j] = directions_oversampled[k];
    BOOST_CHECK(abs(screen_oversampled(i, j).real() - solutions[k].real()) <
                TOLERANCE);
    BOOST_CHECK(abs(screen_oversampled(i, j).imag() - solutions[k].imag()) <
                TOLERANCE);
  }
}
BOOST_AUTO_TEST_SUITE_END()