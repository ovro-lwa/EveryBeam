#ifndef OSKAR_ELEMENTRESPONSE_H
#define OSKAR_ELEMENTRESPONSE_H

#include "../ElementResponse.h"
#include "../Singleton.h"

#include "OSKARDatafile.h"

#include <memory>

namespace everybeam {

//! Implementation of the OSKAR dipole response model
class OSKARElementResponseDipole : public ElementResponse {
 public:
  static std::shared_ptr<OSKARElementResponseDipole> getInstance() {
    return Singleton<OSKARElementResponseDipole>::getInstance();
  }

  virtual void response(
      double freq, double theta, double phi,
      std::complex<double> (&response)[2][2]) const final override;
};

//! Implementation of the OSKAR spherical wave response model
class OSKARElementResponseSphericalWave : public ElementResponse {
 public:
  static std::shared_ptr<OSKARElementResponseSphericalWave> getInstance() {
    return Singleton<OSKARElementResponseSphericalWave>::getInstance();
  }

  OSKARElementResponseSphericalWave();

  virtual void response(
      double freq, double theta, double phi,
      std::complex<double> (&response)[2][2]) const final override;

  virtual void response(
      int element_id, double freq, double theta, double phi,
      std::complex<double> (&response)[2][2]) const final override;

 protected:
  std::string get_path(const char*) const;

  std::unique_ptr<Datafile> m_datafile;
};

}  // namespace everybeam
#endif
