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
/**
 * @brief A BeamFormer contains a number of antennas - be it lower level
 * beamformers or elements - and can return its combined response or array
 * factor.
 *
 */
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
  BeamFormer(const CoordinateSystem& coordinate_system,
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
             const vector3r_t& phase_reference_position,
             std::shared_ptr<FieldResponse> field_response = nullptr)
      : Antenna(coordinate_system, phase_reference_position),
        local_phase_reference_position_(
            TransformToLocalPosition(phase_reference_position_)),
        field_response_(field_response) {}

  BeamFormer(const vector3r_t& phase_reference_position,
             std::shared_ptr<FieldResponse> field_response = nullptr)
      : Antenna(phase_reference_position),
        local_phase_reference_position_(
            TransformToLocalPosition(phase_reference_position_)),
        field_response_(field_response) {}

  std::shared_ptr<Antenna> Clone() const override;

  /**
   * @brief Add an antenna to the antennas_ array.
   *
   * @param antenna
   */
  void AddAntenna(std::shared_ptr<Antenna> antenna) {
    antennas_.push_back(antenna);
    delta_phase_reference_positions_.push_back(
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
  std::shared_ptr<Antenna> ExtractAntenna(size_t antenna_index) const;

  /**
   * @brief Compute the geometric response given the the phase reference
   * directions in the beam former and a direction of interest. In typical use
   * cases, the direction of interest is computed as the (frequency weighted)
   * difference between the pointing direction and the direction of interest,
   * i.e. direction = pointing_freq * pointing_dir - interest_freq *
   *
   * @param phase_reference_positions Phase reference positions.
   * @param direction The direction of interest.
   * @return The geometry response for each position.
   */
  static aocommon::UVector<std::complex<double>> ComputeGeometricResponse(
      const std::vector<vector3r_t>& phase_reference_positions,
      const vector3r_t& direction);

 protected:
  // Compute the BeamFormer response in certain direction of arrival (ITRF, m)
  // and return (Jones) matrix of response
  aocommon::MC2x2 LocalResponse(real_t time, real_t freq,
                                const vector3r_t& direction,
                                const Options& options) const override;

  // Compute the local ArrayFactor, with ArrayFactor a vectorial
  // "representation" of a diagonal Jones matrix
  aocommon::MC2x2Diag LocalArrayFactor(real_t time, real_t freq,
                                       const vector3r_t& direction,
                                       const Options& options) const override;

  const vector3r_t
      local_phase_reference_position_;  // in coordinate system of Antenna

  // List of antennas in BeamFormer
  std::vector<std::shared_ptr<Antenna>> antennas_;
  std::vector<vector3r_t> delta_phase_reference_positions_;

  // (Optional) shared pointer to field response model (e.g. LOBES)
  std::shared_ptr<FieldResponse> field_response_;

 private:
  /**
   * @brief Transform position vector into a local position vector.
   */
  vector3r_t TransformToLocalPosition(const vector3r_t& position);

  /**
   * @brief Compute the beamformer weights based on the difference vector
   * between the pointing direction and the direction of interest. Analogous to
   * \c ComputeGeometricResponse , this difference vector should be computed as:
   * direction = pointing_freq * pointing_dir - interest_freq * interest_dir
   *
   * @param direction Direction of interest (ITRF)
   * @return std::vector<aocommon::MC2x2Diag> Weight matrix per antenna inside
   * the beamformer
   */
  std::vector<aocommon::MC2x2Diag> ComputeWeightedResponses(
      const vector3r_t& direction) const;

  /**
   * @brief Mutex providing a locking mechanism in case field quantities can be
   * precomputed (e.g. the basefunctions in a specific direction for the LOBEs
   * model)
   */
  mutable std::mutex mtx_;
};
}  // namespace everybeam
#endif
