// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "config.h"
#include "../aterms/klfitter.h"
#include <xtensor-fftw/basic.hpp>

#include <vector>
#include <complex>
#include <math.h>

using everybeam::aterms::KlFitter;

#define TOLERANCE 1e-5

BOOST_AUTO_TEST_SUITE(tklfitter)

BOOST_AUTO_TEST_CASE(test1) {
  const int kSubgridSize = 8;
  const int kOrder = 3;

  // Set of directions at exact integer pixels
  std::vector<std::pair<int, int>> directions_int = {{0, 0}, {1, 3}, {4, 6}};

  std::vector<std::pair<float, float>> directions(directions_int.begin(),
                                                  directions_int.end());

  std::vector<size_t> shape = {kSubgridSize, kSubgridSize};
  xt::xarray<float, xt::layout_type::row_major> screen(shape);

  KlFitter kl_fitter(kSubgridSize, kOrder, directions);

  const std::vector<float> solutions = {1.0, 2.0, 3.0};

  kl_fitter.Evaluate(solutions, screen.data());

  for (size_t k = 0; k < directions.size(); ++k) {
    auto [i, j] = directions_int[k];
    BOOST_CHECK(abs(screen(i, j) - solutions[k]) < TOLERANCE);
  }
}

BOOST_AUTO_TEST_SUITE_END()