#include "mwa.h"
#include "../griddedresponse/mwagrid.h"

#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::MWAGrid;
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
    const coords::CoordinateSystem &coordinate_system) {
  std::unique_ptr<GriddedResponse> grid(new MWAGrid(this, coordinate_system));
  return grid;
}
