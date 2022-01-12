// skamidpoint.h: Class for computing the directional telescope
// responses for SKA-MID
//
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_POINTRESPONSE_SKAMIDPOINT_H_
#define EVERYBEAM_POINTRESPONSE_SKAMIDPOINT_H_

#include "dishpoint.h"
#include "../elementresponse.h"

namespace everybeam {
namespace skamidbeam {
class SkaMidResponse;
}

namespace pointresponse {

class SkaMidPoint final : public DishPoint {
 public:
  SkaMidPoint(const telescope::Telescope* telescope_ptr, double time,
              ElementResponseModel element_response_model);

  void Response(BeamMode beam_mode, std::complex<float>* buffer, double ra,
                double dec, double freq, size_t station_idx,
                size_t field_id) override;

 private:
  ElementResponseModel element_response_model_;

  std::unique_ptr<everybeam::skamidbeam::SkaMidResponse> ska_mid_response_;
};
}  // namespace pointresponse
}  // namespace everybeam

#endif