// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "beamformerlofar.h"

#include "common/constants.h"
#include "common/mathutils.h"

#include <cmath>
#include <cassert>

namespace everybeam {
std::vector<std::complex<double>> BeamFormerLofar::ComputeGeometricResponse(
    std::vector<vector3r_t> phase_reference_positions,
    const vector3r_t &direction) const {
  // Initialize and fill result vector by looping over phase_reference_positions
  std::vector<std::complex<double>> result(phase_reference_positions.size());

  for (std::size_t idx = 0; idx < phase_reference_positions.size(); ++idx) {
    // Simplified dot product, since local_phase_reference_position_ = [0, 0, 0]
    const double dl = direction[0] * phase_reference_positions[idx][0] +
                      direction[1] * phase_reference_positions[idx][1] +
                      direction[2] * phase_reference_positions[idx][2];

    // Note that the frequency weighting is already implicit in "direction"
    double phase = -2 * M_PI * dl / common::c;
    result[idx] = {std::cos(phase), std::sin(phase)};
  }
  return result;
}

aocommon::MC2x2Diag BeamFormerLofar::FieldArrayFactor(
    real_t time, real_t freq, const vector3r_t &direction,
    const Options &options, const std::vector<vector3r_t> &antenna_positions,
    const std::vector<std::array<bool, 2>> &antenna_enabled) const {
  assert(antenna_positions.size() == antenna_enabled.size());

  // Weighted subtraction of the directions, with weights given by corresponding
  // freqs. Purpose is to correctly handle the case in which options.freq0 !=
  // freq
  vector3r_t delta_direction =
      options.freq0 * options.station0 - freq * direction;

  // Get geometric response for pointing direction
  std::vector<std::complex<double>> geometric_response =
      ComputeGeometricResponse(antenna_positions, delta_direction);

  double weight_sum[2] = {0.0, 0.0};
  aocommon::MC2x2Diag result(0., 0.);

  for (std::size_t idx = 0; idx < antenna_positions.size(); ++idx) {
    result[0] += geometric_response[idx] * (1.0 * antenna_enabled[idx][0]);
    result[1] += geometric_response[idx] * (1.0 * antenna_enabled[idx][1]);
    weight_sum[0] += (1.0 * antenna_enabled[idx][0]);
    weight_sum[1] += (1.0 * antenna_enabled[idx][1]);
  }

  // Normalize the weight by the number of enabled tiles
  result[0] /= weight_sum[0];
  result[1] /= weight_sum[1];

  return result;
}

aocommon::MC2x2 BeamFormerLofar::LocalResponse(real_t time, real_t freq,
                                               const vector3r_t &direction,
                                               const Options &options) const {
  // Compute the combined array factor
  aocommon::MC2x2Diag array_factor =
      LocalArrayFactor(time, freq, direction, options);

  // NOTE: there are maybe some redundant transformations in element-> response
  aocommon::MC2x2 element_response =
      element_->Response(time, freq, direction, options);
  return array_factor * element_response;
}

}  // namespace everybeam
