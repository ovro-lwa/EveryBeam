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

  void CalculateStation(std::complex<float>* buffer, double ra, double dec,
                        double freq, size_t station_idx,
                        size_t field_id) override;

  void CalculateAllStations(std::complex<float>* buffer, double ra, double dec,
                            double freq, size_t field_id);

 private:
  std::unique_ptr<everybeam::mwabeam::TileBeam2016> tile_beam_;

  mutable std::mutex mtx_;
};
}  // namespace pointresponse
}  // namespace everybeam

#endif  // EVERYBEAM_POINTRESPONSE_MWAPOINT_H_