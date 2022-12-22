// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "config.h"
#include "../aterms/fourierfittingaterm.h"

#include <aocommon/uvector.h>
#include <aocommon/imagecoordinates.h>
#include <vector>
#include <complex>
#include <math.h>

using aocommon::CoordinateSystem;
using everybeam::aterms::FourierFittingATerm;

BOOST_AUTO_TEST_SUITE(tfourierfittingaterm)

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

BOOST_AUTO_TEST_CASE(calculate) {
  std::vector<std::string> station_names = {"CS001HBA0", "CS001HBA1",
                                            "RS509HBA"};

  // double frequency = 122215271.;
  // double ra(-1.44194878), dec(0.85078091);
  double ra(2.2031291147), dec(1.125737367);

  // Properties of grid
  std::size_t width(32), height(32);
  double dl(3.5 * M_PI / 180. / width), dm(3.5 * M_PI / 180. / height),
      shift_l(0.), shift_m(0.);

  CoordinateSystem coord_system = {width, height, ra,      dec,
                                   dl,    dm,     shift_l, shift_m};

  const int support = 3;
  FourierFittingATerm fourierfittingaterm(station_names, coord_system, support);

  // double time = 4987958432.341753;

  aocommon::UVector<std::complex<float>> aterm_buffer(width * height *
                                                      station_names.size() * 4);
  // fourierfittingaterm.Open({"solutions.h5"});
  // fourierfittingaterm.Calculate(aterm_buffer.data(), time, frequency, 0,
  //                               nullptr);
}
BOOST_AUTO_TEST_SUITE_END()
