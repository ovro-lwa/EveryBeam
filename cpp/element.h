// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ELEMENT_H
#define EVERYBEAM_ELEMENT_H

#include <complex>
#include <memory>

#include "antenna.h"
#include "elementresponse.h"
#include "common/types.h"

namespace everybeam {

/**
 * @brief Elementary antenna, for which a response can be computed,
 * but without any substructure like a beamformer
 *
 */
class Element : public Antenna {
 public:
  /**
   * @brief Construct a new Element object
   *
   * @param coordinate_system (antenna) CoordinateSystem
   * @param element_response ElementResponseModel
   * @param id
   */
  Element(const CoordinateSystem& coordinate_system,
          ElementResponse::Ptr element_response, int id)
      : Antenna(coordinate_system),
        id_(id),
        element_response_(element_response) {}

  std::shared_ptr<Antenna> Clone() const override;

  /**
   * @brief Get the Element ID object
   *
   * @return size_t
   */
  size_t GetElementID() const { return id_; }

  /**
   * @brief Convenience function to compute the %Element Response a for given
   * element index
   *
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency of the plane wave (Hz).
   * @param direction Direction of arrival (ITRF, m).
   * @param id Element index.
   * @param options
   * @return aocommon::MC2x2 Jones matrix
   */
  aocommon::MC2x2 ResponseID(real_t time, real_t freq,
                             const vector3r_t& direction, size_t id,
                             const Options& options = {}) {
    // Transform direction and directions in options to local coordinatesystem
    vector3r_t local_direction = TransformToLocalDirection(direction);
    Options local_options;
    local_options.freq0 = options.freq0;
    local_options.station0 = TransformToLocalDirection(options.station0);
    local_options.tile0 = TransformToLocalDirection(options.tile0);
    local_options.rotate = options.rotate;
    local_options.east = TransformToLocalDirection(options.east);
    local_options.north = TransformToLocalDirection(options.north);
    return LocalResponse(time, freq, local_direction, id, local_options);
  }

  /**
   * @brief Compute the local response of the element.
   *
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency of the plane wave (Hz).
   * @param direction Direction of arrival (East-North-Up, m).
   * @param id ID of element
   * @param options
   * @return aocommon::MC2x2
   */
  virtual aocommon::MC2x2 LocalResponse(real_t time, real_t freq,
                                        const vector3r_t& direction, size_t id,
                                        const Options& options) const;

  aocommon::MC2x2Diag ArrayFactor(
      [[maybe_unused]] real_t time, [[maybe_unused]] real_t freq,
      [[maybe_unused]] const vector3r_t& direction,
      [[maybe_unused]] const Options& options) const final override {
    return aocommon::MC2x2Diag::Unity();
  };

 protected:
  aocommon::MC2x2 LocalResponse(real_t time, real_t freq,
                                const vector3r_t& direction,
                                const Options& options) const override {
    return LocalResponse(time, freq, direction, id_, options);
  };

  int id_;
  ElementResponse::Ptr element_response_;
};
}  // namespace everybeam

#endif
