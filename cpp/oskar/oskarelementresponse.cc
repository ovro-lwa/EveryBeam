#include "oskarelementresponse.h"
#include "oskar.h"
#include "config.h"
#include <iostream>

namespace everybeam {

void OSKARElementResponseDipole::Response(
    double freq, double theta, double phi,
    std::complex<double> (&response)[2][2]) const {
  double dipole_length_m = 1;  // TODO
  std::complex<double>* response_ptr = (std::complex<double>*)response;

  double phi_x = phi;
  double phi_y = phi + M_PI_2;
  oskar_evaluate_dipole_pattern_double(1, &theta, &phi_x, freq, dipole_length_m,
                                       response_ptr);
  oskar_evaluate_dipole_pattern_double(1, &theta, &phi_y, freq, dipole_length_m,
                                       response_ptr + 2);
}

OSKARElementResponseSphericalWave::OSKARElementResponseSphericalWave() {
  std::string path = GetPath("oskar.h5");
  datafile_.reset(new Datafile(path));
}

OSKARElementResponseSphericalWave::OSKARElementResponseSphericalWave(
    const std::string& path) {
  datafile_.reset(new Datafile(path));
}

void OSKARElementResponseSphericalWave::Response(
    double freq, double theta, double phi,
    std::complex<double> (&response)[2][2]) const {
  // This ElementResponse model is element specific, so an element_id is
  // required to know for what element the response needs to be evaluated A
  // std::invalid_argument exception is thrown although strictly speaking it are
  // not the given arguments that are invalid, but the Response(...) method with
  // a different signature should have been called.
  throw std::invalid_argument(
      "OSKARElementResponseSphericalWave: missing argument element_id");
}

void OSKARElementResponseSphericalWave::Response(
    int element_id, double freq, double theta, double phi,
    std::complex<double> (&response)[2][2]) const {
  auto dataset = datafile_->Get(freq);
  auto l_max = dataset->GetLMax();

  std::complex<double>* response_ptr = (std::complex<double>*)response;
  std::complex<double>* alpha_ptr = dataset->GetAlphaPtr(element_id);

  double phi_x = phi;
  double phi_y = phi + M_PI_2;
  oskar_evaluate_spherical_wave_sum_double(1, &theta, &phi_x, &phi_y, l_max,
                                           alpha_ptr, response_ptr);
}

std::string OSKARElementResponseSphericalWave::GetPath(
    const char* filename) const {
  std::stringstream ss;
  ss << EVERYBEAM_DATA_DIR << "/";
  ss << filename;
  return ss.str();
}

}  // namespace everybeam
