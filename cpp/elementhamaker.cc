// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "elementhamaker.h"

#include <aocommon/matrix2x2.h>

#include "common/constants.h"
#include "common/mathutils.h"

namespace everybeam {
std::shared_ptr<Antenna> ElementHamaker::Clone() const {
  auto element_clone =
      std::make_shared<ElementHamaker>(coordinate_system_, id_);
  element_clone->enabled_[0] = enabled_[0];
  element_clone->enabled_[1] = enabled_[1];
  return element_clone;
}

aocommon::MC2x2 ElementHamaker::LocalResponse(
    const ElementResponse& element_response, [[maybe_unused]] real_t time,
    real_t freq, const vector3r_t& direction, size_t id,
    const Options& options) const {
  vector2r_t thetaphi = cart2thetaphi(direction);
  thetaphi[1] -= 5.0 * common::pi_4;

  aocommon::MC2x2 result =
      element_response.Response(id, freq, thetaphi[0], thetaphi[1]);

  if (options.rotate) {
    // cross with unit upward pointing vector {0.0, 0.0, 1.0}
    const vector3r_t e_phi = normalize(cross(direction));
    const vector3r_t e_theta = cross(e_phi, direction);
    result *= {dot(e_theta, options.north), dot(e_theta, options.east),
               dot(e_phi, options.north), dot(e_phi, options.east)};
  }
  return result;
}
}  // namespace everybeam
