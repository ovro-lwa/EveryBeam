// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "beamformer.h"

#include "common/constants.h"
#include "common/mathutils.h"

#include <cmath>
#include <cassert>

namespace everybeam {

std::shared_ptr<Antenna> BeamFormer::Clone() const {
  auto beamformer_clone = std::make_shared<BeamFormer>(
      coordinate_system_, phase_reference_position_);

  // antennas_ is a vector of pointers to Antennas, so
  // this creates a shallow copy, in the sense that
  // the antennas are not copied, only the pointers.
  beamformer_clone->antennas_ = antennas_;
  beamformer_clone->delta_phase_reference_positions_ =
      delta_phase_reference_positions_;
  return beamformer_clone;
}

std::shared_ptr<Antenna> BeamFormer::ExtractAntenna(
    size_t antenna_index) const {
  std::shared_ptr<Antenna> antenna = antennas_[antenna_index]->Clone();
  antenna->Transform(coordinate_system_);
  return antenna;
}

vector3r_t BeamFormer::TransformToLocalPosition(const vector3r_t& position) {
  // Get antenna position relative to coordinate system origin
  const vector3r_t dposition{position[0] - coordinate_system_.origin[0],
                             position[1] - coordinate_system_.origin[1],
                             position[2] - coordinate_system_.origin[2]};
  // Return inner product on orthogonal unit vectors of coordinate system
  return {
      dot(coordinate_system_.axes.p, dposition),
      dot(coordinate_system_.axes.q, dposition),
      dot(coordinate_system_.axes.r, dposition),
  };
}

aocommon::UVector<std::complex<double>> BeamFormer::ComputeGeometricResponse(
    const std::vector<vector3r_t>& phase_reference_positions,
    const vector3r_t& direction) {
  constexpr double two_pi_over_c = -2.0 * M_PI / common::c;

  // Allocate and fill result vector by looping over antennas
  aocommon::UVector<std::complex<double>> result(
      phase_reference_positions.size());
  std::vector<double> dl(phase_reference_positions.size());
  std::vector<double> sin_phase(phase_reference_positions.size());
  std::vector<double> cos_phase(phase_reference_positions.size());
#pragma GCC ivdep
  for (size_t i = 0; i < phase_reference_positions.size(); ++i) {
    dl[i] = two_pi_over_c * dot(direction, phase_reference_positions[i]);
  }
// Note that sincos() does not vectorize yet, and
// separate sin() cos() is merged to sincos() by the compiler.
// Hence split the loop into separate sin(), cos() loops.
#pragma GCC ivdep
  for (size_t i = 0; i < phase_reference_positions.size(); ++i) {
    sin_phase[i] = std::sin(dl[i]);
  }
#pragma GCC ivdep
  for (size_t i = 0; i < phase_reference_positions.size(); ++i) {
    cos_phase[i] = std::cos(dl[i]);
  }
#pragma GCC ivdep
  for (size_t i = 0; i < phase_reference_positions.size(); ++i) {
    result[i] = {cos_phase[i], sin_phase[i]};
  }
  return result;
}

std::vector<aocommon::MC2x2Diag> BeamFormer::ComputeWeightedResponses(
    const vector3r_t& pointing) const {
  // Get geometric response for pointing direction
  aocommon::UVector<std::complex<double>> geometric_response =
      ComputeGeometricResponse(delta_phase_reference_positions_, pointing);

  // Initialize and fill result
  double weight_sum[2] = {0.0, 0.0};
  std::vector<aocommon::MC2x2Diag> result(antennas_.size());
  for (size_t idx = 0; idx < antennas_.size(); ++idx) {
    // Get geometric response at index
    const std::complex<double> phasor = geometric_response[idx];
    // Compute the delays in x/y direction
    result[idx] = {phasor * (1.0 * antennas_[idx]->enabled_[0]),
                   phasor * (1.0 * antennas_[idx]->enabled_[1])};
    weight_sum[0] += antennas_[idx]->enabled_[0];
    weight_sum[1] += antennas_[idx]->enabled_[1];
  }

  // Normalize the weight by the number of antennas
  for (auto& entry : result) {
    entry[0] /= weight_sum[0];
    entry[1] /= weight_sum[1];
  }
  return result;
}

aocommon::MC2x2 BeamFormer::LocalResponse(
    const ElementResponse& element_response, real_t time, real_t freq,
    const vector3r_t& direction, const Options& options) const {
  // Weighted subtraction of the pointing direction (0-direction), and the
  // direction of interest. Weights are given by corresponding freqs.
  const vector3r_t delta_direction =
      options.freq0 * options.station0 - freq * direction;

  // Weights based on (weighted) difference vector between
  // pointing direction and direction of interest of beam
  const std::vector<aocommon::MC2x2Diag> weights =
      ComputeWeightedResponses(delta_direction);

  // Copy options into local_options. Needed to propagate
  // the potential change in the rotate boolean downstream.
  Options local_options = options;

  // If fixate_direction_ is true, compute and cache quantities related to the
  // field. This is done for LOBEs beamformers in which all elements inside the
  // beamformer have the same basisfunction for a given direction.
  std::shared_ptr<ElementResponse> local_element_response;
  if (fixate_direction_) {
    local_element_response = element_response.FixateDirection(direction);
    local_options.rotate = false;
  }

  aocommon::MC2x2 result(0.0, 0.0, 0.0, 0.0);
  for (size_t idx = 0; idx < antennas_.size(); ++idx) {
    aocommon::MC2x2 antenna_response = antennas_[idx]->Response(
        local_element_response ? *local_element_response : element_response,
        time, freq, direction, local_options);
    result += weights[idx] * antenna_response;
  }

  // If the Jones matrix needs to be rotated from theta, phi directions
  // to north, east directions, but this has not been done yet, do it here
  if (options.rotate && !local_options.rotate) {
    // cross with unit upward pointing vector {0.0, 0.0, 1.0}
    const vector3r_t e_phi = normalize(cross(direction));
    const vector3r_t e_theta = cross(e_phi, direction);
    result *= {dot(e_theta, options.north), dot(e_theta, options.east),
               dot(e_phi, options.north), dot(e_phi, options.east)};
  }

  return result;
}

aocommon::MC2x2Diag BeamFormer::LocalArrayFactor(real_t time, real_t freq,
                                                 const vector3r_t& direction,
                                                 const Options& options) const {
  // Weighted subtraction of the pointing direction (0-direction), and the
  // direction of interest (direction). Weights are given by corresponding
  // freqs.
  const vector3r_t delta_direction =
      options.freq0 * options.station0 - freq * direction;

  // Weights based on (weighted) difference vector between
  // pointing direction and direction of interest of beam
  const std::vector<aocommon::MC2x2Diag> weights =
      ComputeWeightedResponses(delta_direction);

  aocommon::MC2x2Diag result(0., 0.);
  for (size_t idx = 0; idx < antennas_.size(); ++idx) {
    const aocommon::MC2x2Diag antenna_array_factor =
        antennas_[idx]->ArrayFactor(time, freq, direction, options);
    result += weights[idx] * antenna_array_factor;
  }
  return result;
}

}  // namespace everybeam
