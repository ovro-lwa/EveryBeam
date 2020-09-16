#include "beamformeridenticalantennas.h"

#include "common/constants.h"
#include "common/mathutils.h"

#include <cmath>

namespace everybeam {

Antenna::Ptr BeamFormerIdenticalAntennas::Clone() const {
  auto beamformer_clone = std::make_shared<BeamFormerIdenticalAntennas>(
      coordinate_system_, phase_reference_position_);
  beamformer_clone->antennas_ = antennas_;
  return beamformer_clone;
}

matrix22c_t BeamFormerIdenticalAntennas::LocalResponse(
    real_t time, real_t freq, const vector3r_t &direction,
    const Options &options) const {
  matrix22c_t result;

  auto antenna = antennas_[0];

  matrix22c_t antenna_response =
      antenna->Response(time, freq, direction, options);

  diag22c_t array_factor = LocalArrayFactor(time, freq, direction, options);

  result[0][0] = array_factor[0] * antenna_response[0][0];
  result[0][1] = array_factor[0] * antenna_response[0][1];
  result[1][0] = array_factor[1] * antenna_response[1][0];
  result[1][1] = array_factor[1] * antenna_response[1][1];

  return result;
}

}  // namespace everybeam
