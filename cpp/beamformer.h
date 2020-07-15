#ifndef EVERYBEAM_BEAMFORMER_H
#define EVERYBEAM_BEAMFORMER_H

#include <complex>
#include <vector>

#include "element.h"
#include "common/types.h"

namespace everybeam {
class BeamFormer : public Antenna {
 public:
  typedef std::shared_ptr<BeamFormer> Ptr;

  /**
   * @brief Construct a new BeamFormer object
   *
   */
  BeamFormer() : Antenna() {
    local_phase_reference_position_ =
        TransformToLocalPosition(phase_reference_position_);
  }

  /**
   * @brief Construct a new BeamFormer object given a coordinate system.
   *
   * @param coordinate_system
   */
  BeamFormer(const CoordinateSystem &coordinate_system)
      : Antenna(coordinate_system) {
    local_phase_reference_position_ =
        TransformToLocalPosition(phase_reference_position_);
  }

  /**
   * @brief Construct a new BeamFormer object given a coordinate system and a
   * phase reference position
   *
   * @param coordinate_system
   * @param phase_reference_position
   */
  BeamFormer(CoordinateSystem coordinate_system,
             vector3r_t phase_reference_position)
      : Antenna(coordinate_system, phase_reference_position) {
    local_phase_reference_position_ =
        TransformToLocalPosition(phase_reference_position_);
  }

  BeamFormer(vector3r_t phase_reference_position)
      : Antenna(phase_reference_position) {
    local_phase_reference_position_ =
        TransformToLocalPosition(phase_reference_position_);
  }

  /**
   * @brief Add an antenna to the m_antenna array.
   *
   * @param antenna
   */
  void AddAntenna(Antenna::Ptr antenna) { antennas_.push_back(antenna); }

 private:
  vector3r_t
      local_phase_reference_position_;  // in coordinate system of Antenna

  // Transform position vector into a local position vector
  vector3r_t TransformToLocalPosition(const vector3r_t &position);

  // Compute the BeamFormer response in certain direction of arrival (ITRF, m)
  // and return (Jones) matrix of response
  virtual matrix22c_t LocalResponse(real_t time, real_t freq,
                                    const vector3r_t &direction,
                                    const Options &options) const override;

  // Compute the local ArrayFactor, with ArrayFactor a vectorial
  // "representation" of Jones matrix
  virtual diag22c_t LocalArrayFactor(real_t time, real_t freq,
                                     const vector3r_t &direction,
                                     const Options &options) const override {
    return {1.0, 1.0};
  }

  // Compute the geometric response for all the antennas in the BeamFormer based
  // on the probing frequency and a specified direction (either pointing dir or
  // dir of interest).
  std::vector<std::complex<double>> ComputeGeometricResponse(
      const double freq, const vector3r_t &direction) const;

  // Compute the weights based on the pointing direction of the beam and the
  // beam reference frequence.
  std::vector<std::pair<std::complex<double>, std::complex<double>>>
  ComputeWeights(const vector3r_t &direction, double freq) const;

  // List of antennas in BeamFormer
  std::vector<Antenna::Ptr> antennas_;
};
}  // namespace everybeam
#endif
