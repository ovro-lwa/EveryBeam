// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mwa.h"
#include "../griddedresponse/mwagrid.h"
#include "../pointresponse/mwapoint.h"

#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::MWAGrid;
using everybeam::pointresponse::MWAPoint;
using everybeam::pointresponse::PointResponse;
using everybeam::telescope::MWA;

MWA::MWA(const casacore::MeasurementSet& ms, const Options& options)
    : Telescope(ms, options) {
  if (GetNrStations() == 0) {
    throw std::runtime_error("No antennae in set");
  }

  casacore::MSAntenna antenna(ms.antenna());
  casacore::MPosition::ScalarColumn antenna_pos_col(
      antenna, antenna.columnName(casacore::MSAntennaEnums::POSITION));
  array_position_ = antenna_pos_col(0);

  casacore::Table mwa_tile_pointing =
      ms.keywordSet().asTable("MWA_TILE_POINTING");
  casacore::ArrayColumn<int> delays_col(mwa_tile_pointing, "DELAYS");
  casacore::Array<int> delays_arr = delays_col(0);
  casacore::Array<int>::contiter delays_arr_ptr = delays_arr.cbegin();
  for (auto& delay : delays_) {
    delay = *delays_arr_ptr;
    ++delays_arr_ptr;
  }
}

std::unique_ptr<GriddedResponse> MWA::GetGriddedResponse(
    const coords::CoordinateSystem& coordinate_system) const {
  return std::make_unique<MWAGrid>(this, coordinate_system);
}

std::unique_ptr<PointResponse> MWA::GetPointResponse(double time) const {
  return std::make_unique<MWAPoint>(this, time);
}