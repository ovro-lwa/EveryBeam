// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// hamakerelementresponse.cc:
// Functions to compute the (idealized) response of a LOFAR
// LBA or HBA dual dipole antenna.

#include <stdexcept>

#include "config.h"

#include "hamakerelementresponse.h"
#include "../common/singleton.h"

#include <aocommon/matrix2x2.h>

namespace everybeam {

std::shared_ptr<HamakerElementResponse>
HamakerElementResponse::GetLbaInstance() {
  return common::Singleton<HamakerElementResponseLBA>::GetInstance();
}

std::shared_ptr<HamakerElementResponse> HamakerElementResponse::GetInstance(
    const std::string& name) {
  if (name.find("LBA") != std::string::npos) {
    return common::Singleton<HamakerElementResponseLBA>::GetInstance();
  }
  if (name.find("HBA") != std::string::npos) {
    return common::Singleton<HamakerElementResponseHBA>::GetInstance();
  }
  throw std::invalid_argument(
      "HamakerElementResponse::GetInstance: name should end in either 'LBA' or "
      "'HBA'");
}

std::string HamakerElementResponse::GetPath(const char* filename) const {
  std::stringstream ss;
  ss << EVERYBEAM_DATA_DIR << "/";
  ss << filename;
  return ss.str();
}

aocommon::MC2x2 HamakerElementResponse::Response(double freq, double theta,
                                                 double phi) const {
  aocommon::MC2x2 response = aocommon::MC2x2::Zero();

  // Clip directions below the horizon.
  if (theta >= M_PI_2) {
    return response;
  }

  const double freq_center = coeffs_->GetFreqCenter();
  const double freq_range = coeffs_->GetFreqRange();
  const unsigned int nHarmonics = coeffs_->Get_nHarmonics();
  const unsigned int nPowerTheta = coeffs_->Get_nPowerTheta();
  const unsigned int nPowerFreq = coeffs_->Get_nPowerFreq();

  // The model is parameterized in terms of a normalized frequency in the
  // range [-1, 1]. The appropriate conversion is taken care of below.
  freq = (freq - freq_center) / freq_range;

  // The variables sign and kappa are used to compute the value of kappa
  // mentioned in the description of the beam model [kappa = (-1)^k * (2 * k
  //+ 1)] incrementally.
  int sign = 1, kappa = 1;

  std::pair<std::complex<double>, std::complex<double>> P;
  std::pair<std::complex<double>, std::complex<double>> Pj;
  std::pair<std::complex<double>, std::complex<double>> Pk;
  for (unsigned int k = 0; k < nHarmonics; ++k) {
    // Compute the (diagonal) projection matrix P for the current harmonic.
    // This requires the evaluation of two polynomials in theta and freq (of
    // degree nPowerTheta in theta and nPowerFreq in freq), one for each
    // element of P. The polynomials are evaluated using Horner's rule.

    // Horner's rule requires backward iteration of the coefficients, so
    // start indexing the block of coefficients at the last element

    // Evaluate the highest order term.
    coeffs_->GetCoefficient(k, nPowerTheta - 1, nPowerFreq - 1, P);

    for (unsigned int i = 0; i < nPowerFreq - 1; ++i) {
      coeffs_->GetCoefficient(k, nPowerTheta - 1, nPowerFreq - i - 2, Pk);
      P.first = P.first * freq + Pk.first;
      P.second = P.second * freq + Pk.second;
    }

    // Evaluate the remaining terms.
    for (unsigned int j = 0; j < nPowerTheta - 1; ++j) {
      coeffs_->GetCoefficient(k, nPowerTheta - j - 2, nPowerFreq - 1, Pj);
      for (unsigned int i = 0; i < nPowerFreq - 1; ++i) {
        coeffs_->GetCoefficient(k, nPowerTheta - j - 2, nPowerFreq - i - 2, Pk);
        Pj.first = Pj.first * freq + Pk.first;
        Pj.second = Pj.second * freq + Pk.second;
      }
      P.first = P.first * theta + Pj.first;
      P.second = P.second * theta + Pj.second;
    }

    // Compute the Jones matrix for the current harmonic, by rotating P over
    const double angle = sign * kappa * phi;
    const double caz = std::cos(angle);
    const double saz = std::sin(angle);

    response[0] += caz * P.first;
    response[1] += -saz * P.second;
    response[2] += saz * P.first;
    response[3] += caz * P.second;

    // Update sign and kappa.
    sign = -sign;
    kappa += 2;
  }
  return response;
}

HamakerElementResponseHBA::HamakerElementResponseHBA() {
  std::string path = GetPath("HamakerHBACoeff.h5");
  coeffs_.reset(new HamakerCoefficients(path));
}

HamakerElementResponseLBA::HamakerElementResponseLBA() {
  std::string path = GetPath("HamakerLBACoeff.h5");
  coeffs_.reset(new HamakerCoefficients(path));
}
}  // namespace everybeam
