// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "dish.h"
#include "../griddedresponse/dishgrid.h"
#include "../pointresponse/dishpoint.h"

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using everybeam::griddedresponse::DishGrid;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::pointresponse::DishPoint;
using everybeam::pointresponse::PointResponse;
using everybeam::telescope::Dish;

Dish::Dish(const casacore::MeasurementSet& ms,
           std::unique_ptr<circularsymmetric::Coefficients> coefficients,
           const Options& options)
    : Telescope(ms, options), dish_coefficients_(std::move(coefficients)) {
  casacore::MSField field_table = ms.field();
  casacore::ArrayColumn<double> pointing_dir_col(
      field_table, casacore::MSField::columnName(casacore::MSField::DELAY_DIR));

  for (std::size_t field_id = 0; field_id != field_table.nrow(); ++field_id) {
    casacore::Array<double> pdir = pointing_dir_col(field_id);
    double pdir_ra = *pdir.cbegin();
    double pdir_dec = *(pdir.cbegin() + 1);
    field_pointing_.emplace_back(pdir_ra, pdir_dec);
  }
  SetIsTimeRelevant(false);
}

std::unique_ptr<GriddedResponse> Dish::GetGriddedResponse(
    const aocommon::CoordinateSystem& coordinate_system) const {
  return std::make_unique<DishGrid>(this, coordinate_system);
}

std::unique_ptr<PointResponse> Dish::GetPointResponse(double time) const {
  return std::make_unique<DishPoint>(this, time);
}
