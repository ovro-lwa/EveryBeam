#include "oskar.h"
#include "../common/mathutils.h"
#include "../common/casautils.h"
#include "../msv3readutils.h"

#include <aocommon/banddata.h>
#include <cassert>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using namespace everybeam;
using namespace everybeam::telescope;
using namespace casacore;

OSKAR::OSKAR(MeasurementSet &ms, const ElementResponseModel model,
             const Options &options)
    : Telescope(ms, model, options) {
  ReadAllStations(ms, model);
}

std::unique_ptr<griddedresponse::GriddedResponse> OSKAR::GetGriddedResponse(
    const coords::CoordinateSystem &coordinate_system){
    // Get and return GriddedResponse ptr
    //   std::unique_ptr<griddedresponse::GriddedResponse> grid(
    //       new griddedresponse::LOFARGrid(this, coordinate_system));
    //   // griddedresponse::GriddedResponse grid(LOFARGrid(this,
    //   coordinate_system));
    //   return grid;
};

Station::Ptr OSKAR::ReadStation(const MeasurementSet &ms, std::size_t id,
                                const ElementResponseModel model) const {
  Station::Ptr station = ReadMSv3Station(ms, id, model);
  return station;
}
