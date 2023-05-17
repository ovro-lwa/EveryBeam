// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_HAMAKER_ELEMENTRESPONSE_H_
#define EVERYBEAM_HAMAKER_ELEMENTRESPONSE_H_

#include "../elementresponse.h"
#include "hamakercoeff.h"

#include <memory>

namespace everybeam {

//! Implementation of the Hamaker response model
class [[gnu::visibility("default")]] HamakerElementResponse
    : public ElementResponse {
 public:
  explicit HamakerElementResponse(const std::string& name);

  ElementResponseModel GetModel() const final override {
    return ElementResponseModel::kHamaker;
  }

  aocommon::MC2x2 Response(double freq, double theta, double phi)
      const final override;

  /**
   * @brief Get instance of Hamaker LBA response
   */
  static std::shared_ptr<HamakerElementResponse> GetLbaInstance();

 private:
  // Since coefficients are equal for all LBA and HBA instances, share them.
  static std::weak_ptr<const HamakerCoefficients> cached_lba_coefficients_;
  static std::weak_ptr<const HamakerCoefficients> cached_hba_coefficients_;
  std::shared_ptr<const HamakerCoefficients> coefficients_;
};

}  // namespace everybeam

#endif
