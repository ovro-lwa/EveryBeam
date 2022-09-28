// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>
#include <iostream>
#include <cstdlib>

#include <oskarelementresponse.h>

#include "../../external/npy.hpp"  // to save arrays in numpy format

int main(int argc, char** argv) {
  everybeam::OSKARElementResponseSphericalWave element_response(
      "oskar-comparison.h5");
  double freq = 50e6;

  int N;
  if (argc == 1) {
    N = 256;
  } else {
    N = atoi(argv[1]);
  }

  std::vector<std::complex<double>> result(N * N * 2 * 2);
  aocommon::MC2x2 response;
  for (int i = 0; i < N; ++i) {
    double x = (2.0 * i) / (N - 1) - 1.0;
    for (int j = 0; j < N; ++j) {
      double y = (2.0 * j) / (N - 1) - 1.0;
      double theta = std::asin(sqrt(x * x + y * y));
      double phi = std::atan2(y, x);
      response = element_response.Response(0, freq, theta, phi);
      response.AssignTo(&result[4 * (i * N + j)]);
    }
  }

  const long unsigned leshape[] = {(long unsigned int)N, (long unsigned int)N,
                                   2, 2};
  npy::SaveArrayAsNumpy("response.npy", false, 4, leshape, result);
}
