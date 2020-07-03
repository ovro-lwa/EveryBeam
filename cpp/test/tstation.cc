#include "./../Station.h"

// TODO: make a test out of this
int main() {
  const everybeam::vector3r_t position = {{1.0, 2.0, 3.0}};

  std::string name = "station0_LBA";
  auto model = everybeam::ElementResponseModel::Hamaker;

  // Create station.
  everybeam::Station::Ptr station(
      new everybeam::Station(name, position, model));

  auto element_response = station->get_element_response();

  double freq = 50e6;
  double theta = 0.0;
  double phi = 0.0;
  std::complex<double> response[2][2];

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

      double az = M_PI - phi;
      double el = M_PI_2 - theta;
      element_response->response(0, freq, theta, phi, result_arr[i][j]);
    }
  }
  return 0;
}
