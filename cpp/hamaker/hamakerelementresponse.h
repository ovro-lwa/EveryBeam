// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_HAMAKER_ELEMENTRESPONSE_H_
#define EVERYBEAM_HAMAKER_ELEMENTRESPONSE_H_

#include "../elementresponse.h"
#include "hamakercoeff.h"

#include <memory>

namespace everybeam {

//! Implementation of the Hamaker response model
class HamakerElementResponse : public ElementResponse {
 public:
  ElementResponseModel GetModel() const final override {
    return ElementResponseModel::kHamaker;
  }
  aocommon::MC2x2 Response(double freq, double theta,
                           double phi) const final override;

  /**
   * @brief Get instance of Hamaker LBA response
   */
  static std::shared_ptr<HamakerElementResponse> GetLbaInstance();

  /**
   * @brief Get instance of HamakerElementResponse, infer type (LBA/HBA)
   * from station name.
   *
   * @param name Station Name
   */
  static std::shared_ptr<HamakerElementResponse> GetInstance(
      const std::string& name);

 protected:
  std::string GetPath(const char*) const;

  std::unique_ptr<HamakerCoefficients> coeffs_;
};

class HamakerElementResponseHBA : public HamakerElementResponse {
 public:
  HamakerElementResponseHBA();
};

class HamakerElementResponseLBA : public HamakerElementResponse {
 public:
  HamakerElementResponseLBA();
};

}  // namespace everybeam

#endif
