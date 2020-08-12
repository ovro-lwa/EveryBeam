#include "./../station.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(tstation)

BOOST_AUTO_TEST_CASE(station) {
  const everybeam::vector3r_t position = {{1.0, 2.0, 3.0}};

  std::string name = "station0_LBA";
  auto model = everybeam::ElementResponseModel::kHamaker;

  // Create station.
  everybeam::Station::Ptr station(
      new everybeam::Station(name, position, model));

  auto element_response = station->GetElementResponse();

  double freq = 50e6;

  constexpr int N = 100;
  std::vector<std::complex<double>> result(N * N * 2 * 2);
  typedef std::complex<double> result_arr_t[N][N][2][2];

  result_arr_t &result_arr = *(result_arr_t *)result.data();

  for (int i = 0; i < N; ++i) {
    double x = (2.0 * i) / (N - 1) - 1.0;
    for (int j = 0; j < N; ++j) {
      double y = (2.0 * j) / (N - 1) - 1.0;
      double theta = asin(sqrt(x * x + y * y));
      double phi = atan2(y, x);
      BOOST_REQUIRE_NO_THROW(
          element_response->Response(0, freq, theta, phi, result_arr[i][j]));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
