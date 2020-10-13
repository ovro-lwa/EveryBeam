#include "elementhamaker.h"
#include "common/mathutils.h"
#include "common/constants.h"

namespace everybeam {
Antenna::Ptr ElementHamaker::Clone() const {
  auto element_clone = std::make_shared<ElementHamaker>(coordinate_system_,
                                                        element_response_, id_);
  element_clone->enabled_[0] = enabled_[0];
  element_clone->enabled_[1] = enabled_[1];
  return element_clone;
}

matrix22c_t ElementHamaker::LocalResponse(real_t time, real_t freq,
                                          const vector3r_t &direction,
                                          size_t id,
                                          const Options &options) const {
  vector2r_t thetaphi = cart2thetaphi(direction);
  thetaphi[1] -= 5.0 * common::pi_4;

  matrix22c_t result;
  static_assert(sizeof(std::complex<double>[2][2]) == sizeof(matrix22c_t),
                "matrix22c_t has incorrect size");
  element_response_->Response(
      id, freq, thetaphi[0], thetaphi[1],
      reinterpret_cast<std::complex<double>(&)[2][2]>(result));

  if (options.rotate) {
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
}  // namespace everybeam
