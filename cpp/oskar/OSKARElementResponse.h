#ifndef OSKAR_ELEMENTRESPONSE_H
#define OSKAR_ELEMENTRESPONSE_H

#include "../ElementResponse.h"
#include "../common/Singleton.h"

#include "OSKARDatafile.h"

#include <memory>

namespace everybeam {

//! Implementation of the OSKAR dipole response model
class OSKARElementResponseDipole : public ElementResponse {
 public:
  static std::shared_ptr<OSKARElementResponseDipole> getInstance() {
    return common::Singleton<OSKARElementResponseDipole>::getInstance();
  }

  virtual void response(
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
  static std::shared_ptr<OSKARElementResponseSphericalWave> getInstance() {
    return common::Singleton<OSKARElementResponseSphericalWave>::getInstance();
  }

  /** Constructor loading the default coefficients file */
  OSKARElementResponseSphericalWave();

  /** Constructor loading a custom coefficients file
   *
   * @param path Path to the coefficients file to load
   */
  OSKARElementResponseSphericalWave(const std::string &path);

  virtual void response(
      double freq, double theta, double phi,
      std::complex<double> (&response)[2][2]) const final override;

  virtual void response(
      int element_id, double freq, double theta, double phi,
      std::complex<double> (&response)[2][2]) const final override;

 protected:
  std::string get_path(const char *) const;

  std::unique_ptr<Datafile> m_datafile;
};

}  // namespace everybeam
#endif
