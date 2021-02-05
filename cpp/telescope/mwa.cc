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

MWA::MWA(const casacore::MeasurementSet &ms, const Options &options)
    : Telescope(ms, options) {
  if (GetNrStations() == 0) throw std::runtime_error("No antennae in set");

  casacore::MSAntenna antenna(ms.antenna());
  casacore::MPosition::ScalarColumn antenna_pos_col(
      antenna, antenna.columnName(casacore::MSAntennaEnums::POSITION));
  ms_properties_.array_position = antenna_pos_col(0);

  casacore::Table mwa_tile_pointing =
      ms.keywordSet().asTable("MWA_TILE_POINTING");
  casacore::ArrayColumn<int> delays_col(mwa_tile_pointing, "DELAYS");
  casacore::Array<int> delays_arr = delays_col(0);
  casacore::Array<int>::contiter delays_arr_ptr = delays_arr.cbegin();
  for (int i = 0; i != 16; ++i) ms_properties_.delays[i] = delays_arr_ptr[i];
}

std::unique_ptr<GriddedResponse> MWA::GetGriddedResponse(
    const coords::CoordinateSystem &coordinate_system) const {
  std::unique_ptr<GriddedResponse> grid(new MWAGrid(this, coordinate_system));
  return grid;
}

std::unique_ptr<PointResponse> MWA::GetPointResponse(double time) const {
  // Get and return PointResponse ptr
  std::unique_ptr<PointResponse> point_response(new MWAPoint(this, time));
  return point_response;
}