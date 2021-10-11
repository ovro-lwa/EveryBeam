// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "antenna.h"

#include "common/mathutils.h"

namespace everybeam {

Antenna::Antenna(const CoordinateSystem &coordinate_system,
                 const vector3r_t &phase_reference_position)
    : coordinate_system_(coordinate_system),
      phase_reference_position_(phase_reference_position),
      enabled_{true, true} {}

Antenna::Antenna(const vector3r_t &phase_reference_position)
    : coordinate_system_(
          {phase_reference_position, CoordinateSystem::identity_axes}),
      phase_reference_position_(phase_reference_position),
      enabled_{true, true} {}

constexpr Antenna::CoordinateSystem::Axes
    Antenna::CoordinateSystem::identity_axes;

constexpr vector3r_t Antenna::CoordinateSystem::zero_origin;

constexpr Antenna::CoordinateSystem Antenna::IdentityCoordinateSystem;

vector3r_t Antenna::TransformToLocalDirection(
    const vector3r_t &direction) const {
  return {
      dot(coordinate_system_.axes.p, direction),
      dot(coordinate_system_.axes.q, direction),
      dot(coordinate_system_.axes.r, direction),
  };
}

void Antenna::Transform(const CoordinateSystem &coordinate_system) {
  coordinate_system_.axes.p =
      coordinate_system_.axes.p[0] * coordinate_system.axes.p +
      coordinate_system_.axes.p[1] * coordinate_system.axes.q +
      coordinate_system_.axes.p[2] * coordinate_system.axes.r;

  coordinate_system_.axes.q =
      coordinate_system_.axes.q[0] * coordinate_system.axes.p +
      coordinate_system_.axes.q[1] * coordinate_system.axes.q +
      coordinate_system_.axes.q[2] * coordinate_system.axes.r;

  coordinate_system_.axes.r =
      coordinate_system_.axes.r[0] * coordinate_system.axes.p +
      coordinate_system_.axes.r[1] * coordinate_system.axes.q +
      coordinate_system_.axes.r[2] * coordinate_system.axes.r;

  coordinate_system_.origin =
      coordinate_system.origin +
      coordinate_system_.origin[0] * coordinate_system.axes.p +
      coordinate_system_.origin[1] * coordinate_system.axes.q +
      coordinate_system_.origin[2] * coordinate_system.axes.r;

  phase_reference_position_ =
      coordinate_system.origin +
      phase_reference_position_[0] * coordinate_system.axes.p +
      phase_reference_position_[1] * coordinate_system.axes.q +
      phase_reference_position_[2] * coordinate_system.axes.r;
}

}  // namespace everybeam
