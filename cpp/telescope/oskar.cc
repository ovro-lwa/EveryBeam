#include "oskar.h"
#include "../griddedresponse/oskargrid.h"
#include "../common/mathutils.h"
#include "../common/casautils.h"
#include "../msv3readutils.h"

#include <aocommon/banddata.h>
#include <cassert>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using casacore::MeasurementSet;
using everybeam::Station;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::OSKARGrid;
using everybeam::telescope::OSKAR;

OSKAR::OSKAR(MeasurementSet &ms, const Options &options)
    : Telescope(ms, options) {
  stations_.resize(nstations_);
  ReadAllStations(ms, options_.element_response_model);

  casacore::ScalarMeasColumn<casacore::MDirection> delay_dir_col(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));

  // Populate struct
  ms_properties_ = {.delay_dir = delay_dir_col(0)};
}

std::unique_ptr<GriddedResponse> OSKAR::GetGriddedResponse(
    const coords::CoordinateSystem &coordinate_system) {
  // Get and return GriddedResponse ptr
  std::unique_ptr<GriddedResponse> grid(new OSKARGrid(this, coordinate_system));
  return grid;
}

Station::Ptr OSKAR::ReadStation(const MeasurementSet &ms, std::size_t id,
                                const ElementResponseModel model) const {
  Station::Ptr station = ReadMSv3Station(ms, id, model);
  return station;
}
