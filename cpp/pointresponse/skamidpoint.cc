// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "skamidpoint.h"
#include "../skamidbeam/skamidresponse.h"
#include "../skamidbeam/skamidanalyticalresponse.h"
#include "../telescope/skamid.h"

using everybeam::skamidbeam::SkaMidAnalyticalResponse;
using everybeam::skamidbeam::SkaMidResponse;

namespace everybeam {
namespace pointresponse {

SkaMidPoint::SkaMidPoint(const telescope::Telescope* telescope_ptr, double time,
                         ElementResponseModel element_response_model)
    : DishPoint(telescope_ptr, time),
      element_response_model_(element_response_model) {
  const telescope::SkaMid& ska_mid =
      static_cast<const telescope::SkaMid&>(*telescope_);

  switch (element_response_model_) {
    case ElementResponseModel::kSkaMidAnalytical:
      ska_mid_response_.reset(new SkaMidAnalyticalResponse(
          ska_mid.GetDiameter(), ska_mid.GetBlockage()));
      break;
    default:
      throw std::runtime_error(
          "Element response model not supported for SKA MID");
  }
}

void SkaMidPoint::Response([[maybe_unused]] BeamMode beam_mode,
                           std::complex<float>* buffer, double ra, double dec,
                           double freq, [[maybe_unused]] size_t station_idx,
                           size_t field_id) {
  const telescope::SkaMid& ska_mid_telescope =
      static_cast<const telescope::SkaMid&>(*telescope_);

  double pdir_ra;
  double pdir_dec;
  std::tie(pdir_ra, pdir_dec) = ska_mid_telescope.GetFieldPointing()[field_id];
  ska_mid_response_->Render(buffer, ra, dec, pdir_ra, pdir_dec, freq);
}

}  // namespace pointresponse
}  // namespace everybeam