// pointresponse.h: Base class for computing the directional telescope
// responses.
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_POINTRESPONSE_POINTRESPONSE_H_
#define EVERYBEAM_POINTRESPONSE_POINTRESPONSE_H_

#include <complex>
#include "../telescope/telescope.h"

namespace everybeam {
namespace pointresponse {

/**
 * @brief Virtual base class to compute the point response
 *
 */
class PointResponse {
 public:
  virtual ~PointResponse() {}

  /**
   * @brief Update the (cached) time
   *
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   */
  void UpdateTime(double time) { time_ = time; }

  /**
   * @brief Get beam response for a given station at a prescribed ra, dec
   * position.
   * Method is made virtual, to admit more efficient overrides in the near
   * future.
   *
   * @param buffer Buffer with a size of 4 complex floats to receive the beam
   * response
   * @param ra Right ascension (rad)
   * @param dec Declination (rad)
   * @param freq Frequency (Hz)
   * @param station_idx Station index
   * @param field_id
   */
  virtual void CalculateStation(std::complex<float>* buffer, double ra,
                                double dec, double freq, size_t station_idx,
                                size_t field_id) = 0;

  /**
   * @brief Same as PointResponseStation, but now iterate over all stations in
   * MS.
   *
   * @param buffer Buffer with a size of 4 * nr_stations complex floats to
   * receive the beam response
   * @param ra Right ascension (rad)
   * @param dec Declination (rad)
   * @param freq Frequency (Hz)
   * @param field_id
   */
  virtual void CalculateAllStations(std::complex<float>* buffer, double ra,
                                    double dec, double freq, size_t field_id) {
    for (size_t i = 0; i < telescope_->GetNrStations(); ++i) {
      CalculateStation(buffer, ra, dec, freq, i, field_id);
      buffer += 4;
    }
  };

  std::size_t GetAllStationsBufferSize() const {
    return telescope_->GetNrStations() * 4u;
  }

 protected:
  /**
   * @brief Construct a new Point Response object
   * NOTE: constructor will become protected once this class forms the
   * a virtual base class for all specializations
   *
   * @param telescope_ptr Const pointer to telescope object
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   */
  PointResponse(const telescope::Telescope* telescope_ptr, double time)
      : telescope_(telescope_ptr), time_(time){};

  const telescope::Telescope* telescope_;
  double time_;
};
}  // namespace pointresponse
}  // namespace everybeam

#endif  // EVERYBEAM_POINTRESPONSE_POINTRESPONSE_H_