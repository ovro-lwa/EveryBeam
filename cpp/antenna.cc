#include "antenna.h"

#include "common/mathutils.h"

namespace everybeam {

vector3r_t Antenna::TransformToLocalDirection(
    const vector3r_t &direction) const {
  vector3r_t local_direction{
      dot(coordinate_system_.axes.p, direction),
      dot(coordinate_system_.axes.q, direction),
      dot(coordinate_system_.axes.r, direction),
  };

  return local_direction;
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
