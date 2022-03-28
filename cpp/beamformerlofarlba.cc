// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "beamformerlofarlba.h"

namespace everybeam {

std::shared_ptr<Antenna> BeamFormerLofarLBA::Clone() const {
  auto beamformer_clone = std::make_shared<BeamFormerLofarLBA>(
      coordinate_system_, phase_reference_position_);

  // NOTE: this is an incomplete clone, only creating a deep-copy of the
  // element. In fact, it also hides an upcast from an ElementHamaker into
  // an Element object.
  // The sole and single purpose of Clone() is to be used in
  // Station::SetAntenna!
  Element element_copy = *element_;
  beamformer_clone->SetElement(std::make_shared<Element>(element_copy));
  return beamformer_clone;
}

aocommon::MC2x2Diag BeamFormerLofarLBA::LocalArrayFactor(
    real_t time, real_t freq, const vector3r_t& direction,
    const Options& options) const {
  // Compute the array factor of the field
  return FieldArrayFactor(time, freq, direction, options, element_positions_,
                          element_enabled_);
}
}  // namespace everybeam