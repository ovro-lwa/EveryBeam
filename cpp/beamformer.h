// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_BEAMFORMER_H
#define EVERYBEAM_BEAMFORMER_H

#include <complex>
#include <vector>
#include <mutex>

#include <aocommon/uvector.h>

#include "element.h"
#include "common/types.h"
#include "common/mathutils.h"
#include "fieldresponse.h"

namespace everybeam {
class BeamFormer : public Antenna {
 public:
  typedef std::shared_ptr<BeamFormer> Ptr;

  /**
   * @brief Construct a new BeamFormer object
   *
   */
  BeamFormer(std::shared_ptr<FieldResponse> field_response = nullptr)
      : Antenna(),
        local_phase_reference_position_(
            TransformToLocalPosition(phase_reference_position_)),
        field_response_(field_response) {}

  /**
   * @brief Construct a new BeamFormer object given a coordinate system.
   *
   * @param coordinate_system
   */
  BeamFormer(const CoordinateSystem &coordinate_system,
             std::shared_ptr<FieldResponse> field_response = nullptr)
      : Antenna(coordinate_system),
        local_phase_reference_position_(
            TransformToLocalPosition(phase_reference_position_)),
        field_response_(field_response) {}

  /**
   * @brief Construct a new BeamFormer object given a coordinate system and a
   * phase reference position
   *
   * @param coordinate_system
   * @param phase_reference_position
   */
  BeamFormer(CoordinateSystem coordinate_system,
             const vector3r_t &phase_reference_position,
             std::shared_ptr<FieldResponse> field_response = nullptr)
      : Antenna(coordinate_system, phase_reference_position),
        local_phase_reference_position_(
            TransformToLocalPosition(phase_reference_position_)),
        field_response_(field_response) {}

  BeamFormer(const vector3r_t &phase_reference_position,
             std::shared_ptr<FieldResponse> field_response = nullptr)
      : Antenna(phase_reference_position),
        local_phase_reference_position_(
            TransformToLocalPosition(phase_reference_position_)),
        field_response_(field_response) {}

  Antenna::Ptr Clone() const override;

  /**
   * @brief Add an antenna to the antennas_ array.
   *
   * @param antenna
   */
  void AddAntenna(Antenna::Ptr antenna) {
    antennas_.push_back(antenna);
    delta_phase_reference_position_.push_back(
        antennas_.back()->phase_reference_position_ -
        local_phase_reference_position_);
  }

  /**
   * @brief Extracts an antenna from the beamformer
   *
   * @param antenna_index index of antenna to extact
   * @returns pointer to a copy of antenna with index antenna_index
   *
   * The antenna is extracted such that it can be used stand-alone,
   * independent of the beamformer. The coordinate system of the extracted
   * antenna is transformed from internal representation to external
   * representation by application of the beamformer coordinate system to
   * the antenna coordinate system.
   *
   * The returned antenna can be either an Element or a BeamFormer.
   *
   * The beamformer itself remains unchanged.
   */
  Antenna::Ptr ExtractAntenna(size_t antenna_index) const;

 protected:
  const vector3r_t
      local_phase_reference_position_;  // in coordinate system of Antenna

  // Transform position vector into a local position vector
  vector3r_t TransformToLocalPosition(const vector3r_t &position);

  // Compute the BeamFormer response in certain direction of arrival (ITRF, m)
  // and return (Jones) matrix of response
  aocommon::MC2x2 LocalResponse(real_t time, real_t freq,
                                const vector3r_t &direction,
                                const Options &options) const override;

  // Compute the local ArrayFactor, with ArrayFactor a vectorial
  // "representation" of Jones matrix
  aocommon::MC2x2Diag LocalArrayFactor(real_t time, real_t freq,
                                       const vector3r_t &direction,
                                       const Options &options) const override;

  // Compute the geometric response for all the antennas in the BeamFormer based
  // on the difference vector between the pointing direction and the direction
  // of interest. This difference vector should be computed as: direction =
  // pointing_freq * pointing_dir - interest_freq * interest_dir
  aocommon::UVector<std::complex<double>> ComputeGeometricResponse(
      const vector3r_t &direction) const;

  // Compute the weights based on the difference vector between the pointing
  // direction and the direction of interest. Analogous to
  // ComputeGeometricResponse, this difference vector should be computed as:
  // direction = pointing_freq * pointing_dir - interest_freq * interest_dir
  std::vector<aocommon::MC2x2Diag> ComputeWeightedResponses(
      const vector3r_t &direction) const;

  // List of antennas in BeamFormer
  std::vector<Antenna::Ptr> antennas_;
  std::vector<vector3r_t> delta_phase_reference_position_;

  // (Optional) shared pointer to field response model (e.g. LOBES)
  std::shared_ptr<FieldResponse> field_response_;

 private:
  mutable std::mutex mtx_;
};
}  // namespace everybeam
#endif
