// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "./../station.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(tstation)

BOOST_AUTO_TEST_CASE(station) {
  const everybeam::vector3r_t position = {{1.0, 2.0, 3.0}};

  std::string name = "station0_LBA";
  everybeam::Options options;
  options.element_response_model = everybeam::ElementResponseModel::kHamaker;

  // Create station.
  std::shared_ptr<everybeam::Station> station =
      std::make_shared<everybeam::Station>(name, position, options);

  const everybeam::ElementResponse& element_response =
      *station->GetElementResponse();

  double freq = 50e6;

  constexpr int N = 100;
  for (int i = 0; i < N; ++i) {
    double x = (2.0 * i) / (N - 1) - 1.0;
    for (int j = 0; j < N; ++j) {
      double y = (2.0 * j) / (N - 1) - 1.0;
      double theta = asin(sqrt(x * x + y * y));
      double phi = atan2(y, x);
      BOOST_REQUIRE_NO_THROW(element_response.Response(0, freq, theta, phi));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
