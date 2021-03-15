// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef OSKAR_ELEMENTRESPONSE_H
#define OSKAR_ELEMENTRESPONSE_H

#include "../elementresponse.h"
#include "../common/singleton.h"

#include "oskardatafile.h"

#include <memory>

namespace everybeam {

//! Implementation of the OSKAR dipole response model
class OSKARElementResponseDipole : public ElementResponse {
 public:
  static std::shared_ptr<OSKARElementResponseDipole> GetInstance() {
    return common::Singleton<OSKARElementResponseDipole>::GetInstance();
  }

  ElementResponseModel GetModel() const final override { return kOSKARDipole; }

  virtual void Response(
      double freq, double theta, double phi,
      std::complex<double> (&response)[2][2]) const final override;
};

//! Implementation of the OSKAR spherical wave response model
class OSKARElementResponseSphericalWave : public ElementResponse {
 public:
  /**
   * A constructor-like static method to instantiate the class
   *
   * returns a globally shared instance of the class that is instantiated
   * in the first call
   */
  static std::shared_ptr<OSKARElementResponseSphericalWave> GetInstance() {
    return common::Singleton<OSKARElementResponseSphericalWave>::GetInstance();
  }

  /** Constructor loading the default coefficients file */
  OSKARElementResponseSphericalWave();

  /** Constructor loading a custom coefficients file
   *
   * @param path Path to the coefficients file to load
   */
  OSKARElementResponseSphericalWave(const std::string &path);

  ElementResponseModel GetModel() const final override {
    return kOSKARSphericalWave;
  }

  virtual void Response(
      double freq, double theta, double phi,
      std::complex<double> (&response)[2][2]) const final override;

  virtual void Response(
      int element_id, double freq, double theta, double phi,
      std::complex<double> (&response)[2][2]) const final override;

 protected:
  std::string GetPath(const char *) const;

  std::unique_ptr<Datafile> datafile_;
};

}  // namespace everybeam
#endif
