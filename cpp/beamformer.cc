#include "beamformer.h"

#include "common/constants.h"
#include "common/mathutils.h"

#include <cmath>

namespace everybeam {

Antenna::Ptr BeamFormer::Clone() const {
  auto beamformer_clone = std::make_shared<BeamFormer>(
      coordinate_system_, phase_reference_position_);

  // antennas_ is a vector of pointers to Antennas, so
  // this creates a shallow copy, in the sense that
  // the antennas are not copied, only the pointers.
  beamformer_clone->antennas_ = antennas_;

  return beamformer_clone;
}

Antenna::Ptr BeamFormer::ExtractAntenna(size_t antenna_index) const {
  Antenna::Ptr antenna = antennas_[antenna_index]->Clone();
  antenna->Transform(coordinate_system_);
  return antenna;
}

vector3r_t BeamFormer::TransformToLocalPosition(const vector3r_t &position) {
  // Get antenna position relative to coordinate system origin
  vector3r_t dposition{position[0] - coordinate_system_.origin[0],
                       position[1] - coordinate_system_.origin[1],
                       position[2] - coordinate_system_.origin[2]};
  // Inner product on orthogonal unit vectors of coordinate system
  vector3r_t local_position{
      dot(coordinate_system_.axes.p, dposition),
      dot(coordinate_system_.axes.q, dposition),
      dot(coordinate_system_.axes.r, dposition),
  };

  return local_position;
}

std::vector<std::complex<double>> BeamFormer::ComputeGeometricResponse(
    const vector3r_t &direction) const {
  // Initialize and fill result vector by looping over antennas
  std::vector<std::complex<double>> result(antennas_.size());
  for (std::size_t idx = 0; idx < antennas_.size(); ++idx) {
    const double dl =
        direction[0] * (antennas_[idx]->phase_reference_position_[0] -
                        local_phase_reference_position_[0]) +
        direction[1] * (antennas_[idx]->phase_reference_position_[1] -
                        local_phase_reference_position_[1]) +
        direction[2] * (antennas_[idx]->phase_reference_position_[2] -
                        local_phase_reference_position_[2]);
    // Note that the frequency is (and should be!) implicit in dl!
    double phase = -2 * M_PI * dl / common::c;
    result[idx] = {std::cos(phase), std::sin(phase)};
  }
  return result;
}

std::vector<std::pair<std::complex<double>, std::complex<double>>>
BeamFormer::ComputeWeightedResponses(const vector3r_t &pointing) const {
  // Get geometric response for pointing direction
  std::vector<std::complex<double>> geometric_response =
      ComputeGeometricResponse(pointing);

  // Initialize and fill result
  double weight_sum[2] = {0.0, 0.0};
  std::vector<std::pair<std::complex<double>, std::complex<double>>> result(
      geometric_response.size());
  for (std::size_t idx = 0; idx < antennas_.size(); ++idx) {
    // Get geometric response at index
    std::complex<double> phasor = geometric_response[idx];
    // Compute the delays in x/y direction
    result[idx] = {phasor * (1.0 * antennas_[idx]->enabled_[0]),
                   phasor * (1.0 * antennas_[idx]->enabled_[1])};
    weight_sum[0] += (1.0 * antennas_[idx]->enabled_[0]);
    weight_sum[1] += (1.0 * antennas_[idx]->enabled_[1]);
  }

  // Normalize the weight by the number of antennas
  for (std::size_t idx = 0; idx < antennas_.size(); ++idx) {
    result[idx].first /= weight_sum[0];
    result[idx].second /= weight_sum[1];
  }

  return result;
}

matrix22c_t BeamFormer::LocalResponse(real_t time, real_t freq,
                                      const vector3r_t &direction,
                                      const Options &options) const {
  std::unique_lock<std::mutex> lock(mtx_, std::defer_lock);

  // Weighted subtraction of the pointing direction (0-direction), and the
  // direction of interest. Weights are given by corresponding freqs.
  vector3r_t delta_direction =
      options.freq0 * options.station0 - freq * direction;

  // Weights based on (weighted) difference vector between
  // pointing direction and direction of interest of beam
  std::vector<std::pair<std::complex<double>, std::complex<double>>> weights =
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
    vector2r_t thetaphi = cart2thetaphi(direction);
    field_response_->SetFieldQuantities(thetaphi[0], thetaphi[1]);
    local_options.rotate = false;
  }

  matrix22c_t result = {{{0}}};
  for (std::size_t antenna_idx = 0; antenna_idx < antennas_.size();
       ++antenna_idx) {
    Antenna::Ptr antenna = antennas_[antenna_idx];
    std::pair<std::complex<double>, std::complex<double>> antenna_weight =
        weights[antenna_idx];

    matrix22c_t antenna_response =
        antenna->Response(time, freq, direction, local_options);

    result[0][0] += antenna_weight.first * antenna_response[0][0];
    result[0][1] += antenna_weight.first * antenna_response[0][1];
    result[1][0] += antenna_weight.second * antenna_response[1][0];
    result[1][1] += antenna_weight.second * antenna_response[1][1];
  }

  // If the Jones matrix needs to be rotated from theta, phi directions
  // to north, east directions, but this has not been done yet, do it here
  if (options.rotate && !local_options.rotate) {
    vector3r_t up = {0.0, 0.0, 1.0};
    vector3r_t e_phi = normalize(cross(up, direction));
    vector3r_t e_theta = cross(e_phi, direction);
    matrix22r_t rotation;
    rotation[0] = {dot(e_theta, options.north), dot(e_theta, options.east)};
    rotation[1] = {dot(e_phi, options.north), dot(e_phi, options.east)};
    result = result * rotation;
  }

  return result;
}

diag22c_t BeamFormer::LocalArrayFactor(real_t time, real_t freq,
                                       const vector3r_t &direction,
                                       const Options &options) const {
  // Weighted subtraction of the pointing direction (0-direction), and the
  // direction of interest (direction). Weights are given by corresponding
  // freqs.
  vector3r_t delta_direction =
      options.freq0 * options.station0 - freq * direction;

  // Weights based on (weighted) difference vector between
  // pointing direction and direction of interest of beam
  std::vector<std::pair<std::complex<double>, std::complex<double>>> weights =
      ComputeWeightedResponses(delta_direction);

  diag22c_t result = {0};
  for (std::size_t antenna_idx = 0; antenna_idx < antennas_.size();
       ++antenna_idx) {
    Antenna::Ptr antenna = antennas_[antenna_idx];
    std::pair<std::complex<double>, std::complex<double>> antenna_weight =
        weights[antenna_idx];

    result[0] += antenna_weight.first;
    result[1] += antenna_weight.second;
  }
  return result;
}

}  // namespace everybeam
