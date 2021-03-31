// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "beamformeridenticalantennas.h"

#include "common/constants.h"
#include "common/mathutils.h"

#include <cmath>

namespace everybeam {

Antenna::Ptr BeamFormerIdenticalAntennas::Clone() const {
  auto beamformer_clone = std::make_shared<BeamFormerIdenticalAntennas>(
      coordinate_system_, phase_reference_position_);
  beamformer_clone->antennas_ = antennas_;
  return beamformer_clone;
}

aocommon::MC2x2 BeamFormerIdenticalAntennas::LocalResponse(
    real_t time, real_t freq, const vector3r_t &direction,
    const Options &options) const {
  auto antenna = antennas_[0];

  aocommon::MC2x2 antenna_response =
      antenna->Response(time, freq, direction, options);
  aocommon::MC2x2Diag array_factor =
      LocalArrayFactor(time, freq, direction, options);
  return array_factor * antenna_response;
}

}  // namespace everybeam
