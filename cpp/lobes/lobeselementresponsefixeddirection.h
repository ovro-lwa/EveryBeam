// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_LOBES_ELEMENTRESPONSE_FIXED_DIRECTION_H
#define EVERYBEAM_LOBES_ELEMENTRESPONSE_FIXED_DIRECTION_H

#include "lobeselementresponse.h"

namespace everybeam {

/**
 * LOBES element response model with a fixed direction.
 *
 * Fixing the direction allows reusing the base functions in different
 * Response() calls.
 */
class LobesElementResponseFixedDirection : public ElementResponse {
 public:
  /**
   * @brief Construct a new LobesElementResponseFixedDirection object that wraps
   *        an existing LOBESElementResponse object.
   *
   * @param element_response The element response object that should be wrapped.
   * @param base_functions The fixed base functions.
   */
  explicit LobesElementResponseFixedDirection(
      std::shared_ptr<const LOBESElementResponse> element_response,
      LOBESElementResponse::BaseFunctions base_functions)
      : element_response_(std::move(element_response)),
        base_functions_(std::move(base_functions)) {}

  ElementResponseModel GetModel() const final override {
    return ElementResponseModel::kLOBES;
  }

  /**
   * @brief Stub override of the Response method, an element id
   * is needed to compute the element response.
   */
  aocommon::MC2x2 Response([[maybe_unused]] double freq,
                           [[maybe_unused]] double theta,
                           [[maybe_unused]] double phi) const final override {
    throw std::invalid_argument(
        "LobesElementResponseFixedDirection::Response needs an element_id");
  };

  /**
   * @brief Virtual implementation of Response method
   *
   * @param element_id ID of element
   * @param freq Frequency of the plane wave (Hz).
   * @param theta Angle wrt. z-axis (rad) (NOTE: parameter is just a stub to
   * match override)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad) (NOTE: parameter is
   * just a stub to match override)
   * @return The Jones matrix as a 2x2 array.
   */
  aocommon::MC2x2 Response(int element_id, double frequency,
                           [[maybe_unused]] double theta,
                           [[maybe_unused]] double phi) const final override {
    // Clip directions below the horizon.
    if (theta >= M_PI_2) {
      return aocommon::MC2x2::Zero();
    }

    return element_response_->Response(base_functions_, element_id, frequency);
  }

 private:
  std::shared_ptr<const LOBESElementResponse> element_response_;
  LOBESElementResponse::BaseFunctions base_functions_;
};
}  // namespace everybeam

#endif
