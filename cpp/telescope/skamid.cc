// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "skamid.h"
#include "../pointresponse/skamidpoint.h"
#include "../griddedresponse/skamidgrid.h"
#include "../skamidbeam/skamidresponse.h"

namespace {
static constexpr double kDiameter = 15.0;
static constexpr double kBlockage = 0.0;
}  // namespace

namespace everybeam {
namespace telescope {

SkaMid::SkaMid(const casacore::MeasurementSet& ms, const Options& options)
    : Telescope(ms, options),
      element_response_model_(options.element_response_model) {
  casacore::MSField field_table = ms.field();
  casacore::ArrayColumn<double> pointing_dir_col(
      field_table, casacore::MSField::columnName(casacore::MSField::DELAY_DIR));

  for (size_t field_id = 0; field_id != field_table.nrow(); ++field_id) {
    casacore::Array<double> pdir = pointing_dir_col(field_id);
    const double pdir_ra = *pdir.cbegin();
    const double pdir_dec = *(pdir.cbegin() + 1);
    field_pointing_.emplace_back(pdir_ra, pdir_dec);
  }

  if (element_response_model_ == ElementResponseModel::kSkaMidAnalytical) {
    SetIsTimeRelevant(false);
  } else {
    SetIsTimeRelevant(true);
  }
}

std::unique_ptr<griddedresponse::GriddedResponse> SkaMid::GetGriddedResponse(
    const aocommon::CoordinateSystem& coordinate_system) const {
  return std::make_unique<griddedresponse::SkaMidGrid>(this, coordinate_system,
                                                       element_response_model_);
}

std::unique_ptr<pointresponse::PointResponse> SkaMid::GetPointResponse(
    double time) const {
  return std::make_unique<pointresponse::SkaMidPoint>(this, time,
                                                      element_response_model_);
}

double SkaMid::GetDiameter() const { return kDiameter; }

double SkaMid::GetBlockage() const { return kBlockage; }

}  // namespace telescope
}  // namespace everybeam
