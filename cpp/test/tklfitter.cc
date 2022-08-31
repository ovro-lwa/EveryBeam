// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../aterms/klfitter.h"

#include <complex>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <xtensor/xarray.hpp>

#include "config.h"

using everybeam::aterms::KlFitter;

BOOST_AUTO_TEST_SUITE(klfitter)

BOOST_AUTO_TEST_CASE(evaluate) {
  constexpr int kSubgridSize = 8;
  constexpr int kOrder = 3;
  constexpr double kTolerance = 1.0e-3;

  // Set of directions at exact integer pixels
  const std::vector<std::pair<int, int>> directions_int = {
      {0, 0}, {1, 3}, {4, 6}};

  const std::vector<std::pair<float, float>> directions(directions_int.begin(),
                                                        directions_int.end());

  const std::vector<size_t> shape = {kSubgridSize, kSubgridSize};
  xt::xarray<float, xt::layout_type::row_major> screen(shape);

  KlFitter kl_fitter(kSubgridSize, kOrder, directions);

  const std::vector<float> solutions = {1.0, 2.0, 3.0};

  kl_fitter.Evaluate(solutions, screen.data());

  for (int i = 0; i < kSubgridSize; ++i) {
    for (int j = 0; j < kSubgridSize; ++j) {
      // Perform some generic checks.
      BOOST_CHECK_GT(screen(i, j), 0.9);
      BOOST_CHECK_LT(screen(i, j), 4.5);
      if (j > 0) {
        BOOST_CHECK_GT(screen(i, j), screen(i, j - 1));
      }

      // At the directions, the screen values should match the solution.
      if (std::make_pair(i, j) == directions_int[0]) {
        BOOST_CHECK_CLOSE(screen(i, j), solutions[0], kTolerance);
      } else if (std::make_pair(i, j) == directions_int[1]) {
        BOOST_CHECK_CLOSE(screen(i, j), solutions[1], kTolerance);
      } else if (std::make_pair(i, j) == directions_int[2]) {
        BOOST_CHECK_CLOSE(screen(i, j), solutions[2], kTolerance);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()