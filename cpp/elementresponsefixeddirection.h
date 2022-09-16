// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ELEMENTRESPONSE_FIXED_DIRECTION_H
#define EVERYBEAM_ELEMENTRESPONSE_FIXED_DIRECTION_H

#include "elementresponse.h"
#include "common/mathutils.h"

namespace everybeam {

/**
 * Element response model with a fixed direction.
 */
class ElementResponseFixedDirection : public ElementResponse {
 public:
  /**
   * @brief Construct a new ElementResponseFixedDirection object that wraps
   *        an existing ElementResponse object.
   *
   * @param element_response The element response object that should be wrapped.
   * @param theta The fixed theta value for all Response() calls.
   * @param phi The fixed phi value for all Response() calls.
   */
  explicit ElementResponseFixedDirection(
      std::shared_ptr<const ElementResponse> element_response, double theta,
      double phi)
      : element_response_(std::move(element_response)),
        theta_(theta),
        phi_(phi) {}

  ElementResponseModel GetModel() const final override {
    return element_response_->GetModel();
  }

  /**
   * @brief Override of the Response method that uses the fixed theta and phi
   * values instead of the supplied arguments.
   */
  aocommon::MC2x2 Response(double freq, [[maybe_unused]] double theta,
                           [[maybe_unused]] double phi) const final override {
    return element_response_->Response(freq, theta_, phi_);
  };

  /**
   * @brief Override of the Response method that uses the fixed theta and phi
   * values instead of the supplied arguments.
   */
  aocommon::MC2x2 Response(int element_id, double frequency,
                           [[maybe_unused]] double theta,
                           [[maybe_unused]] double phi) const final override {
    return element_response_->Response(element_id, frequency, theta_, phi_);
  }

  /**
   * @brief Avoids putting a ElementResponseFixedDirection
   * on top of the current ElementResponseFixedDirection object.
   */
  std::shared_ptr<ElementResponse> FixateDirection(
      const vector3r_t& direction) const final override {
    const vector2r_t thetaphi = cart2thetaphi(direction);
    return std::make_shared<ElementResponseFixedDirection>(
        element_response_, thetaphi[0], thetaphi[1]);
  }

 private:
  std::shared_ptr<const ElementResponse> element_response_;
  double theta_;
  double phi_;
};
}  // namespace everybeam

#endif
