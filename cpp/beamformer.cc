// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "beamformer.h"

#include "common/constants.h"
#include "common/mathutils.h"

#include <cmath>
#include <cassert>

namespace everybeam {

Antenna::Ptr BeamFormer::Clone() const {
  auto beamformer_clone = std::make_shared<BeamFormer>(
      coordinate_system_, phase_reference_position_);

  // antennas_ is a vector of pointers to Antennas, so
  // this creates a shallow copy, in the sense that
  // the antennas are not copied, only the pointers.
  beamformer_clone->antennas_ = antennas_;
  beamformer_clone->delta_phase_reference_position_ =
      delta_phase_reference_position_;
  return beamformer_clone;
}

Antenna::Ptr BeamFormer::ExtractAntenna(size_t antenna_index) const {
  Antenna::Ptr antenna = antennas_[antenna_index]->Clone();
  antenna->Transform(coordinate_system_);
  return antenna;
}

vector3r_t BeamFormer::TransformToLocalPosition(const vector3r_t &position) {
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
    const vector3r_t &direction) const {
  constexpr double two_pi_over_c = -2.0 * M_PI / common::c;

  // Allocate and fill result vector by looping over antennas
  assert(antennas_.size() == delta_phase_reference_position_.size());
  aocommon::UVector<std::complex<double>> result(antennas_.size());
  for (size_t idx = 0; idx < delta_phase_reference_position_.size(); ++idx) {
    const double dl = dot(direction, delta_phase_reference_position_[idx]);
    // Note that the frequency is (and should be!) implicit in dl!
    // We could save a multiplication here, by pre-multiplying the
    // delta_phase_reference_position by two_pi_over_c
    const double phase = two_pi_over_c * dl;
    result[idx] = {std::cos(phase), std::sin(phase)};
  }
  return result;
}

std::vector<aocommon::MC2x2Diag> BeamFormer::ComputeWeightedResponses(
    const vector3r_t &pointing) const {
  // Get geometric response for pointing direction
  aocommon::UVector<std::complex<double>> geometric_response =
      ComputeGeometricResponse(pointing);

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
  for (auto &entry : result) {
    entry[0] /= weight_sum[0];
    entry[1] /= weight_sum[1];
  }
  return result;
}

aocommon::MC2x2 BeamFormer::LocalResponse(real_t time, real_t freq,
                                          const vector3r_t &direction,
                                          const Options &options) const {
  std::unique_lock<std::mutex> lock(mtx_, std::defer_lock);

  // Weighted subtraction of the pointing direction (0-direction), and the
  // direction of interest. Weights are given by corresponding freqs.
  const vector3r_t delta_direction =
      options.freq0 * options.station0 - freq * direction;

  // Weights based on (weighted) difference vector between
  // pointing direction and direction of interest of beam
  const std::vector<aocommon::MC2x2Diag> weights =
      ComputeWeightedResponses(delta_direction);

  // Copy options into local_options. Needed to propagate
  // the potential change in the rotate boolean downstream
  Options local_options = options;
  // If field_response_ not nullptr, set/precompute quantities
  // related to the field
  if (nullptr != field_response_.get()) {
    // Lock the associated mutex, thus avoiding that the LOBESElementResponse
    // basefunctions_ are overwritten before response is computed
    lock.lock();
    const vector2r_t thetaphi = cart2thetaphi(direction);
    field_response_->SetFieldQuantities(thetaphi[0], thetaphi[1]);
    local_options.rotate = false;
  }

  aocommon::MC2x2 result(0.0, 0.0, 0.0, 0.0);
  for (size_t idx = 0; idx < antennas_.size(); ++idx) {
    aocommon::MC2x2 antenna_response =
        antennas_[idx]->Response(time, freq, direction, local_options);
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

  // Wipe out basefunctions cache
  if (nullptr != field_response_.get()) {
    field_response_->ClearFieldQuantities();
  }
  return result;
}

aocommon::MC2x2Diag BeamFormer::LocalArrayFactor(real_t time, real_t freq,
                                                 const vector3r_t &direction,
                                                 const Options &options) const {
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
