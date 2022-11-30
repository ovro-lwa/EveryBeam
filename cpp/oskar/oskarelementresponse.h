// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef OSKAR_ELEMENTRESPONSE_H
#define OSKAR_ELEMENTRESPONSE_H

#include "../elementresponse.h"

#include <memory>

namespace everybeam {

// Use a forward declaration instead of including internal libeverybeam-oskar
// headers. The internal headers use OSKAR headers, which are not available
// outside the library.
class Datafile;

//! Implementation of the OSKAR dipole response model
class OSKARElementResponseDipole : public ElementResponse {
 public:
  ElementResponseModel GetModel() const final override {
    return ElementResponseModel::kOSKARDipole;
  }

  aocommon::MC2x2 Response(double freq, double theta,
                           double phi) const final override;
};

//! Implementation of the OSKAR spherical wave response model
class OSKARElementResponseSphericalWave : public ElementResponse {
 public:
  /** Constructor loading the default coefficients file */
  OSKARElementResponseSphericalWave();

  /**
   * Constructor loading a custom coefficients file
   * @param filename Filename of HDF5 file with coefficients.
   */
  OSKARElementResponseSphericalWave(const std::string& filename);

  ElementResponseModel GetModel() const final override {
    return ElementResponseModel::kOSKARSphericalWave;
  }

  aocommon::MC2x2 Response(double freq, double theta,
                           double phi) const final override;

  aocommon::MC2x2 Response(int element_id, double freq, double theta,
                           double phi) const final override;

 private:
  // This weak pointer allows reusing the default coefficients.
  static std::weak_ptr<Datafile> cached_datafile_;
  std::shared_ptr<Datafile> datafile_;
};

}  // namespace everybeam
#endif
