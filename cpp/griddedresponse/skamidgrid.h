// skamidgrid.h: Class for computing the circular symmetric (gridded) response.
//
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_SKAMIDGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_SKAMIDGRID_H_

#include "dishgrid.h"
#include "../elementresponse.h"

namespace everybeam {

namespace skamidbeam {
class SkaMidResponse;
}

namespace griddedresponse {

class SkaMidGrid final : public DishGrid {
 public:
  SkaMidGrid(const telescope::Telescope* telescope_ptr,
             const aocommon::CoordinateSystem coordinate_system,
             ElementResponseModel element_response_model);

  void Response(BeamMode beam_mode, std::complex<float>* buffer, double time,
                double frequency, size_t station_idx, size_t field_id) override;

 private:
  ElementResponseModel element_response_model_;

  std::unique_ptr<everybeam::skamidbeam::SkaMidResponse> ska_mid_response_;
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_SKAMIDGRID_H_
