// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_BEAMFORMERLOFAR_H
#define EVERYBEAM_BEAMFORMERLOFAR_H

#include <complex>
#include <vector>
#include <memory>

#include "beamformer.h"
#include "elementhamaker.h"
#include "common/types.h"
#include "antenna.h"

namespace everybeam {
/**
 * @brief Abstract class that serves as the base class for the optimized
 * implementations of the beam formers for the LOFAR HBA and LBA telescope in
 * combination with the Hamaker element response model.
 *
 */
class BeamFormerLofar : public Antenna {
 public:
  /**
   * @brief Construct a new BeamFormerLofar object
   *
   */
  BeamFormerLofar() : Antenna() {}

  /**
   * @brief Construct a new BeamFormerLofar object given a coordinate system.
   *
   * @param coordinate_system
   */
  BeamFormerLofar(const CoordinateSystem& coordinate_system)
      : Antenna(coordinate_system) {}

  /**
   * @brief Construct a new BeamFormerLofar object given a coordinate system
   * and a phase reference position
   *
   * @param coordinate_system
   * @param phase_reference_position
   */
  BeamFormerLofar(CoordinateSystem coordinate_system,
                  vector3r_t phase_reference_position)
      : Antenna(coordinate_system, phase_reference_position) {}

  BeamFormerLofar(vector3r_t phase_reference_position)
      : Antenna(phase_reference_position) {}

  /**
   * @brief Pure virtual method of clone.
   *
   * @return std::shared_ptr<Antenna>
   */
  std::shared_ptr<Antenna> Clone() const override = 0;

  /**
   * @brief Set the (unique) Element object for the BeamFormerLofar object.
   *
   * @param element
   */
  void SetElement(std::shared_ptr<Element> element) { element_ = element; }

  /**
   * @brief Add element position to the element_positions_ array
   *
   * @param position
   */
  void AddElementPosition(const vector3r_t& position) {
    element_positions_.push_back(position);
  }

  /**
   * @brief Returns the element_ object of the BeamFormerLofarHBA
   *
   * @return std::shared_ptr<Element>
   */
  std::shared_ptr<Element> GetElement() const { return element_; };

 protected:
  // Pure virtual implementation of array factor at station level
  aocommon::MC2x2Diag LocalArrayFactor(
      real_t time, real_t freq, const vector3r_t& direction,
      const Options& options) const override = 0;

  // Array factor at Field level. antenna_positions_ and antenna_enabled_
  // either represent the tiles in case of LOFAR HBA or the elements in case
  // of LOFAR LBA
  aocommon::MC2x2Diag FieldArrayFactor(
      real_t time, real_t freq, const vector3r_t& direction,
      const Options& options, const std::vector<vector3r_t>& antenna_positions,
      const std::vector<std::array<bool, 2>>& antenna_enabled) const;

  aocommon::MC2x2 LocalResponse(const ElementResponse& element_response,
                                real_t time, real_t freq,
                                const vector3r_t& direction,
                                const Options& options) const final override;

  // Each BeamFormerLofar stores one unique element, and a vector of unique
  // positions. Usually, the Element will be of type ElementHamaker.
  std::shared_ptr<Element> element_;
  std::vector<vector3r_t> element_positions_;
};
}  // namespace everybeam
#endif