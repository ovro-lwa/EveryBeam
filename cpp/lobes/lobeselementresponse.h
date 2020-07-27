#ifndef LOBES_ELEMENTRESPONSE_H
#define LOBES_ELEMENTRESPONSE_H

#include "../elementresponse.h"

#include <memory>

namespace everybeam {

//! Implementation of the Lobes response model
class LOBESElementResponse : public ElementResponse {
 public:
  LOBESElementResponse(std::string name);

  virtual void Response(
      double freq, double theta, double phi,
      std::complex<double> (&response)[2][2]) const final override;

  static std::shared_ptr<LOBESElementResponse> GetInstance(std::string name);
};

}  // namespace everybeam

#endif