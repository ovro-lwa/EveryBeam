#ifndef EVERYBEAM_BEAMFORMERIDENTICALANTENNAS_H
#define EVERYBEAM_BEAMFORMERIDENTICALANTENNAS_H

#include "beamformer.h"

namespace everybeam {
class BeamFormerIdenticalAntennas : public BeamFormer {
 public:
  /**
   * @brief Construct a new BeamFormerIdenticalAntennas object
   *
   */
  BeamFormerIdenticalAntennas() : BeamFormer() {}

  /**
   * @brief Construct a new BeamFormerIdenticalAntennas object given a
   * coordinate system.
   *
   * @param coordinate_system
   */
  BeamFormerIdenticalAntennas(const CoordinateSystem &coordinate_system)
      : BeamFormer(coordinate_system) {}

  /**
   * @brief Construct a new BeamFormer object given a coordinate system and a
   * phase reference position
   *
   * @param coordinate_system
   * @param phase_reference_position
   */
  BeamFormerIdenticalAntennas(CoordinateSystem coordinate_system,
                              vector3r_t phase_reference_position)
      : BeamFormer(coordinate_system, phase_reference_position) {}

  BeamFormerIdenticalAntennas(vector3r_t phase_reference_position)
      : BeamFormer(phase_reference_position) {}

 private:
  // Compute the BeamFormer response in certain direction of arrival (ITRF, m)
  // and return (Jones) matrix of response
  virtual matrix22c_t LocalResponse(real_t time, real_t freq,
                                    const vector3r_t &direction,
                                    const Options &options) const override;
};
}  // namespace everybeam
#endif
