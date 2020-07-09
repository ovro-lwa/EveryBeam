#include "antenna.h"

#include "common/math_utils.h"

namespace everybeam {
vector3r_t Antenna::TransformToLocalDirection(const vector3r_t &direction) {
  vector3r_t local_direction{
      dot(m_coordinate_system.axes.p, direction),
      dot(m_coordinate_system.axes.q, direction),
      dot(m_coordinate_system.axes.r, direction),
  };

  return local_direction;
}
}  // namespace everybeam