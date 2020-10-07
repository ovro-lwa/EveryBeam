#include "beamformerlofarhba.h"

#include "common/constants.h"
#include "common/mathutils.h"

#include <cmath>

namespace everybeam {

Antenna::Ptr BeamFormerLofarHBA::Clone() const {
  auto beamformer_clone = std::make_shared<BeamFormerLofarHBA>(
      BeamFormerLofarHBA(coordinate_system_, phase_reference_position_));

  // NOTE: this is an incomplete clone, only creating a deep-copy of the
  // element. In fact, it also hides an upcast from an ElementHamaker into
  // an Element object.
  // The sole and single purpose of Clone() is to be used in
  // Station::SetAntenna!
  Element element_copy = *element_;
  beamformer_clone->SetElement(std::make_shared<Element>(element_copy));
  return beamformer_clone;
}

std::vector<std::complex<double>> BeamFormerLofarHBA::ComputeGeometricResponse(
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

diag22c_t BeamFormerLofarHBA::LocalArrayFactor(real_t time, real_t freq,
                                               const vector3r_t &direction,
                                               const Options &options) const {
  diag22c_t result = {0};

  // Compute the array factor of the field
  diag22c_t array_factor_field =
      FieldArrayFactor(time, freq, direction, options);

  // Compute the array factor of a tile
  std::complex<double> array_factor_tile =
      TileArrayFactor(time, freq, direction, options);

  result[0] = array_factor_tile * array_factor_field[0];
  result[1] = array_factor_tile * array_factor_field[1];

  return result;
}

diag22c_t BeamFormerLofarHBA::FieldArrayFactor(real_t time, real_t freq,
                                               const vector3r_t &direction,
                                               const Options &options) const {
  // Weighted subtraction of the directions, with weights given by corresponding
  // freqs. Purpose is to correctly handle the case in which options.freq0 !=
  // freq
  vector3r_t delta_direction =
      options.freq0 * options.station0 - freq * direction;

  // Get geometric response for pointing direction
  std::vector<std::complex<double>> geometric_response =
      ComputeGeometricResponse(tile_positions_, delta_direction);

  double weight_sum[2] = {0.0, 0.0};
  diag22c_t result;

  for (std::size_t idx = 0; idx < tile_positions_.size(); ++idx) {
    result[0] += geometric_response[idx] * (1.0 * tile_enabled_[idx][0]);
    result[1] += geometric_response[idx] * (1.0 * tile_enabled_[idx][1]);
    weight_sum[0] += (1.0 * tile_enabled_[idx][0]);
    weight_sum[1] += (1.0 * tile_enabled_[idx][1]);
  }

  // Normalize the weight by the number of enabled tiles
  result[0] /= weight_sum[0];
  result[1] /= weight_sum[1];

  return result;
}

std::complex<double> BeamFormerLofarHBA::TileArrayFactor(
    real_t time, real_t freq, const vector3r_t &direction,
    const Options &options) const {
  // Weighted subtraction of the directions, with weights given by corresponding
  // freqs. Purpose is to correctly handle the case in which options.freq0 !=
  // freq
  vector3r_t delta_direction = options.freq0 * options.tile0 - freq * direction;

  // Get geometric response for the difference vector stored in "pointing"
  std::vector<std::complex<double>> geometric_response =
      ComputeGeometricResponse(element_positions_, delta_direction);

  // Initialize and fill result
  std::complex<double> result = 0;
  for (std::size_t idx = 0; idx < element_positions_.size(); ++idx) {
    result += geometric_response[idx];
  }

  // Normalize the result by the number of tiles
  double weight = element_positions_.size();
  result /= weight;

  return result;
}

matrix22c_t BeamFormerLofarHBA::LocalResponse(real_t time, real_t freq,
                                              const vector3r_t &direction,
                                              const Options &options) const {
  matrix22c_t result = {0};

  // Compute the combined array factor
  diag22c_t array_factor = LocalArrayFactor(time, freq, direction, options);

  // NOTE: there are maybe some redundant transformations in element-> response
  matrix22c_t element_response =
      element_->Response(time, freq, direction, options);

  result[0][0] = (array_factor[0] * element_response[0][0]);
  result[0][1] = (array_factor[0] * element_response[0][1]);
  result[1][0] = (array_factor[1] * element_response[1][0]);
  result[1][1] = (array_factor[1] * element_response[1][1]);

  return result;
}
}  // namespace everybeam
