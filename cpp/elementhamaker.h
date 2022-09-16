// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ELEMENT_HAMAKER_H
#define EVERYBEAM_ELEMENT_HAMAKER_H

#include <complex>
#include <memory>

#include "antenna.h"
#include "element.h"
#include "elementresponse.h"
#include "common/types.h"

namespace everybeam {

/**
 * @brief Elementary antenna, optimized for LOFAR Hamaker model. Derived from
 * the Element class.
 *
 */
class ElementHamaker final : public Element {
 public:
  /**
   * @brief Construct a new Element object
   *
   * @param coordinate_system (antenna) CoordinateSystem
   * @param element_response ElementResponseModel
   * @param id
   */
  ElementHamaker(const CoordinateSystem& coordinate_system, int id)
      : Element(coordinate_system, id) {}

  std::shared_ptr<Antenna> Clone() const override;

  /**
   * @brief This override avoids a number of redundant coordinate
   * transformations compared to the parent implementation.
   */
  aocommon::MC2x2 Response(const ElementResponse& element_response, real_t time,
                           real_t freq, const vector3r_t& direction,
                           const Options& options) const override {
    // The only transform that is needed is hard-coded in LocalResponse
    return LocalResponse(element_response, time, freq, direction, options);
  }

  aocommon::MC2x2 LocalResponse(const ElementResponse& element_response,
                                real_t time, real_t freq,
                                const vector3r_t& direction, size_t id,
                                const Options& options) const override;

 private:
  aocommon::MC2x2 LocalResponse(const ElementResponse& element_response,
                                real_t time, real_t freq,
                                const vector3r_t& direction,
                                const Options& options) const override {
    return LocalResponse(element_response, time, freq, direction, id_, options);
  };
};
}  // namespace everybeam

#endif
