// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "beamformerlofarhba.h"
#include "common/mathutils.h"

namespace everybeam {

Antenna::Ptr BeamFormerLofarHBA::Clone() const {
  auto beamformer_clone = std::make_shared<BeamFormerLofarHBA>(
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

aocommon::MC2x2Diag BeamFormerLofarHBA::LocalArrayFactor(
    real_t time, real_t freq, const vector3r_t &direction,
    const Options &options) const {
  // Compute the array factor of the field
  aocommon::MC2x2Diag array_factor_field = FieldArrayFactor(
      time, freq, direction, options, tile_positions_, tile_enabled_);

  // Compute the array factor of a tile
  const std::complex<double> array_factor_tile =
      TileArrayFactor(time, freq, direction, options);

  return array_factor_field * array_factor_tile;
}

std::complex<double> BeamFormerLofarHBA::TileArrayFactor(
    [[maybe_unused]] real_t time, real_t freq, const vector3r_t &direction,
    const Options &options) const {
  // Weighted subtraction of the directions, with weights given by corresponding
  // freqs. Purpose is to correctly handle the case in which options.freq0 !=
  // freq
  vector3r_t delta_direction = options.freq0 * options.tile0 - freq * direction;

  // Get geometric response for the difference vector stored in "pointing"
  const aocommon::UVector<std::complex<double>> geometric_response =
      ComputeGeometricResponse(element_positions_, delta_direction);

  // Initialize and fill result
  std::complex<double> result = 0;
  for (const auto &gr : geometric_response) {
    result += gr;
  }

  // Normalize the result by the number of tiles
  const double weight = element_positions_.size();
  result /= weight;

  return result;
}
}  // namespace everybeam
