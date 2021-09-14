// dishpoint.h: Class for computing a circular symmetric beam response at given
// point
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_POINTRESPONSE_DISHPOINT_H_
#define EVERYBEAM_POINTRESPONSE_DISHPOINT_H_

#include "pointresponse.h"

namespace everybeam {
namespace pointresponse {

/**
 * @brief Class for computing the directional response of dish telescopes,
 * e.g. VLA, ATCA.
 *
 */
class DishPoint final : public PointResponse {
 public:
  DishPoint(const telescope::Telescope* telescope_ptr, double time)
      : PointResponse(telescope_ptr, time){};

  void FullBeam(std::complex<float>* buffer, double ra, double dec, double freq,
                size_t station_idx, size_t field_id) override;

  void FullBeamAllStations(std::complex<float>* buffer, double ra, double dec,
                           double freq, size_t field_id) override;
};
}  // namespace pointresponse
}  // namespace everybeam
#endif  // EVERYBEAM_POINTRESPONSE_DISHPOINT_H_