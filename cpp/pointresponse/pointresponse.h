// pointresponse.h: Base class for computing the directional telescope
// responses.
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_POINTRESPONSE_POINTRESPONSE_H_
#define EVERYBEAM_POINTRESPONSE_POINTRESPONSE_H_

#include <complex>
#include <mutex>

#include "../common/types.h"
#include "../telescope/telescope.h"

#include <aocommon/matrix2x2diag.h>
#include <aocommon/matrix2x2.h>

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
  void UpdateTime(double time) {
    // Second condition enables "backwards marching" in time
    if (time - time_ > update_interval_ || time_ - time > 0.0) {
      time_ = time;
      has_time_update_ = true;
    } else {
      has_time_update_ = false;
    }
  }

  /**
   * @brief Set interval for updating the time. Can be used for caching
   * ITRF direction vectors.
   *
   * @param update_interval Update interval (s)
   */
  void SetUpdateInterval(double update_interval) {
    update_interval_ = update_interval;
    has_time_update_ = true;
  }

  /**
   * @brief Check whether cached time settings have changed
   */
  bool HasTimeUpdate() const { return has_time_update_; }

  /**
   * @brief See FullBeam.
   */
  [[deprecated("Use FullBeam() instead.")]] void CalculateStation(
      std::complex<float>* response_matrix, double ra, double dec, double freq,
      size_t station_idx, size_t field_id) {
    FullBeam(response_matrix, ra, dec, freq, station_idx, field_id);
  };

  /**
   * @brief See FullBeamAllStations.
   */
  [[deprecated("Use FullBeamAllStations() instead.")]] virtual void
  CalculateAllStations(std::complex<float>* response_matrices, double ra,
                       double dec, double freq, size_t field_id) {
    FullBeamAllStations(response_matrices, ra, dec, freq, field_id);
  };

  /**
   * @brief Get beam response for a given station at a prescribed ra, dec
   * position.
   *
   * @param response_matrix Buffer with a size of 4 complex floats to receive
   * the beam response
   * @param ra Right ascension (rad)
   * @param dec Declination (rad)
   * @param freq Frequency (Hz)
   * @param station_id Station index, corresponding to measurement set antenna
   * index.
   * @param field_id Field index as used in the measurement set
   */
  virtual void FullBeam(std::complex<float>* response_matrix, double ra,
                        double dec, double freq, size_t station_id,
                        size_t field_id) = 0;

  /**
   * @brief Get the full beam response for a station, given a pointing direction
   * in ITRF coordinates
   *
   * @param station_idx Station index
   * @param freq Frequency (Hz)
   * @param direction Direction in ITRF
   * @param mutex Optional mutex. When provided, the caller keeps control over
   * thread-safety. If not provided, the internal mutex will be used and the
   * caller is assumed to be thread-safe.
   * @return aocommon::MC2x2
   */
  virtual aocommon::MC2x2 FullBeam(size_t station_idx, double freq,
                                   const vector3r_t& direction,
                                   std::mutex* mutex = nullptr) {
    throw std::runtime_error("Not yet implemented");
  };

  /**
   * @brief Same as FullBeam, but now iterate over all stations in
   * MS.
   *
   * @param response_matrices Buffer with a size of 4 * nr_stations complex
   * floats to receive the beam response
   * @param ra Right ascension (rad)
   * @param dec Declination (rad)
   * @param freq Frequency (Hz)
   * @param field_id Field index as used in the measurement set
   */
  virtual void FullBeamAllStations(std::complex<float>* response_matrices,
                                   double ra, double dec, double freq,
                                   size_t field_id) {
    for (size_t i = 0; i < telescope_->GetNrStations(); ++i) {
      FullBeam(response_matrices, ra, dec, freq, i, field_id);
      response_matrices += 4;
    }
  };

  /**
   * @brief Get the array factor for a station, given a pointing direction
   * in ITRF coordinates
   *
   * @param station_idx Station index
   * @param freq Frequency (Hz)
   * @param direction Direction in ITRF
   * @param mutex Optional mutex. When provided, the caller keeps control over
   * thread-safety. If not provided, the internal mutex will be used and the
   * caller is assumed to be thread-safe.
   * @return aocommon::MC2x2Diag
   */
  virtual aocommon::MC2x2Diag ArrayFactor(size_t station_idx, double freq,
                                          const vector3r_t& direction,
                                          std::mutex* mutex = nullptr) {
    throw std::runtime_error("Not yet implemented");
  };

  virtual aocommon::MC2x2 ElementResponse(size_t station_idx, double freq,
                                          const vector3r_t& direction) const {
    throw std::runtime_error("Not yet implemented");
  }

  std::size_t GetAllStationsBufferSize() const {
    return telescope_->GetNrStations() * 4u;
  }

 protected:
  /**
   * @brief Construct a new Point Response object
   *
   * @param telescope_ptr Const pointer to telescope object
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   */
  PointResponse(const telescope::Telescope* telescope_ptr, double time)
      : telescope_(telescope_ptr),
        time_(time),
        update_interval_(0),
        has_time_update_(true){};

  const telescope::Telescope* telescope_;
  double time_;
  double update_interval_;
  bool has_time_update_;
};
}  // namespace pointresponse
}  // namespace everybeam

#endif  // EVERYBEAM_POINTRESPONSE_POINTRESPONSE_H_
