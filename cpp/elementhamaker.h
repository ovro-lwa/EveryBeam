// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
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
  ElementHamaker(const CoordinateSystem &coordinate_system,
                 ElementResponse::Ptr element_response, int id)
      : Element(coordinate_system, element_response, id) {}

  std::shared_ptr<Antenna> Clone() const override;

  aocommon::MC2x2 Response(real_t time, real_t freq,
                           const vector3r_t &direction,
                           const Options &options) const override {
    // The only transform that is needed is hard-coded in LocalResponse
    return LocalResponse(time, freq, direction, options);
  }

  aocommon::MC2x2 LocalResponse(real_t time, real_t freq,
                                const vector3r_t &direction, size_t id,
                                const Options &options) const override;

 private:
  aocommon::MC2x2 LocalResponse(real_t time, real_t freq,
                                const vector3r_t &direction,
                                const Options &options) const override {
    return LocalResponse(time, freq, direction, id_, options);
  };
};
}  // namespace everybeam

#endif
