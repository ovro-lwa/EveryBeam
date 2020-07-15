#include "antenna.h"

#include "common/math_utils.h"

namespace everybeam {
vector3r_t Antenna::TransformToLocalDirection(const vector3r_t &direction) {
  vector3r_t local_direction{
      dot(coordinate_system_.axes.p, direction),
      dot(coordinate_system_.axes.q, direction),
      dot(coordinate_system_.axes.r, direction),
  };

  return local_direction;
}
}  // namespace everybeam