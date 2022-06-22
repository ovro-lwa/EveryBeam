// mwapoint.h: Class for computing the MWA beam response at given
// point
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_POINTRESPONSE_MWAPOINT_H_
#define EVERYBEAM_POINTRESPONSE_MWAPOINT_H_

#include "pointresponse.h"
#include "../mwabeam/tilebeam2016.h"

#include <mutex>

namespace everybeam {
namespace pointresponse {

class MWAPoint final : public PointResponse {
 public:
  MWAPoint(const telescope::Telescope* telescope_ptr, double time)
      : PointResponse(telescope_ptr, time){};

  /**
   * @brief Get beam response for a given station at a prescribed ra, dec
   * position.
   * NOTE: function complies with the standard
   * threading rules, but does not guarantee thread-safety itself for efficiency
   * reasons. The caller is responsible to ensure this.
   *
   * @param buffer Buffer with a size of 4 complex floats to receive the beam
   * response
   * @param ra Right ascension (rad)
   * @param dec Declination (rad)
   * @param freq Frequency (Hz)
   * @param station_idx Station index
   * @param field_id
   */
  void Response(BeamMode beam_mode, std::complex<float>* buffer, double ra,
                double dec, double freq, size_t station_idx,
                size_t field_id) override;

  void ResponseAllStations(BeamMode beam_mode, std::complex<float>* buffer,
                           double ra, double dec, double freq,
                           size_t field_id) override;

 private:
  void SetJ200Vectors();
  std::unique_ptr<everybeam::mwabeam::TileBeam2016> tile_beam_;

  casacore::MDirection::Ref j2000_ref_;
  casacore::MDirection::Convert j2000_to_hadecref_;
  casacore::MDirection::Convert j2000_to_azelgeoref_;

  double arr_latitude_;
  std::mutex mutex_;
};
}  // namespace pointresponse
}  // namespace everybeam

#endif  // EVERYBEAM_POINTRESPONSE_MWAPOINT_H_
