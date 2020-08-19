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

OSKAR::OSKAR(MeasurementSet &ms, const Options &options)
    : Telescope(ms, options) {
  stations_.resize(nstations_);
  ReadAllStations(ms, options_.element_response_model);
}

std::unique_ptr<griddedresponse::GriddedResponse> OSKAR::GetGriddedResponse(
    const coords::CoordinateSystem &coordinate_system) {
  throw std::runtime_error(
      "GetGriddedResponse() is not implemented for OSKAR Telescope");

  // TODO: return an OSKARGrid here, in a similar way to the commented out code
  // below

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
