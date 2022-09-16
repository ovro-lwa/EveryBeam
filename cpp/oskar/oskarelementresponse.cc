// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oskarelementresponse.h"

#include <iostream>

#include "config.h"

#include <oskar_beam_utils.h>

using oskar::beam_utils::oskar_evaluate_dipole_pattern_double;
using oskar::beam_utils::oskar_evaluate_spherical_wave_sum_double;

namespace everybeam {

aocommon::MC2x2 OSKARElementResponseDipole::Response(double freq, double theta,
                                                     double phi) const {
  aocommon::MC2x2 response = aocommon::MC2x2::Zero();
  double dipole_length_m = 1;  // TODO

  double phi_x = phi;
  double phi_y = phi + M_PI_2;
  oskar_evaluate_dipole_pattern_double(1, &theta, &phi_x, freq, dipole_length_m,
                                       response.Data());
  oskar_evaluate_dipole_pattern_double(1, &theta, &phi_y, freq, dipole_length_m,
                                       response.Data() + 2);
  return response;
}

OSKARElementResponseSphericalWave::OSKARElementResponseSphericalWave() {
  std::string path = GetPath("oskar.h5");
  datafile_.reset(new Datafile(path));
}

OSKARElementResponseSphericalWave::OSKARElementResponseSphericalWave(
    const std::string& path) {
  datafile_.reset(new Datafile(path));
}

aocommon::MC2x2 OSKARElementResponseSphericalWave::Response(
    [[maybe_unused]] double freq, [[maybe_unused]] double theta,
    [[maybe_unused]] double phi) const {
  // This ElementResponse model is element specific, so an element_id is
  // required to know for what element the response needs to be evaluated A
  // std::invalid_argument exception is thrown although strictly speaking it are
  // not the given arguments that are invalid, but the Response(...) method with
  // a different signature should have been called.
  throw std::invalid_argument(
      "OSKARElementResponseSphericalWave: missing argument element_id");
}

aocommon::MC2x2 OSKARElementResponseSphericalWave::Response(int element_id,
                                                            double freq,
                                                            double theta,
                                                            double phi) const {
  aocommon::MC2x2 response = aocommon::MC2x2::Zero();

  const Dataset& dataset = datafile_->Get(freq);
  const size_t l_max = dataset.GetLMax();

  const std::complex<double>* alpha_ptr = dataset.GetAlphaPtr(element_id);

  double phi_x = phi;
  double phi_y = phi;

  // TODO: phi_x and phi_y can have different values if there is only one set
  // of coefficients that is is used for both dipoles.
  // In that case it is assumed the Y dipole rotated 90deg with respect
  // to the X dipole, so then phi_y = phi+ M_PI_2.
  // That case needs to be detected when the coefficients are read,
  // and here phi_y needs to be set accordingly.

  oskar_evaluate_spherical_wave_sum_double(theta, phi_x, phi_y, l_max,
                                           alpha_ptr, response.Data());
  return response;
}

std::string OSKARElementResponseSphericalWave::GetPath(
    const char* filename) const {
  std::stringstream ss;
  ss << EVERYBEAM_DATA_DIR << "/";
  ss << filename;
  return ss.str();
}

}  // namespace everybeam
