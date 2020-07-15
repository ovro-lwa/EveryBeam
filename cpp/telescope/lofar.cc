#include "lofar.h"
#include "./../gridded_response/lofargrid.h"
#include "./../common/math_utils.h"
#include "./../common/casa_utils.h"
#include "../LofarMetaDataUtil.h"

#include <aocommon/banddata.h>
#include <cassert>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using namespace everybeam;
using namespace everybeam::telescope;
using namespace casacore;

LOFAR::LOFAR(MeasurementSet &ms, const ElementResponseModel model,
             const Options &options)
    : Telescope(ms, model, options) {
  ReadAllStations(ms, model);

  // Populate MeasurementSet properties struct
  aocommon::BandData band(ms.spectralWindow());
  MSAntenna antenna(ms.antenna());
  MPosition::ScalarColumn antenna_pos_col(
      antenna, antenna.columnName(MSAntennaEnums::POSITION));
  MEpoch::ScalarColumn time_column(ms,
                                   ms.columnName(casacore::MSMainEnums::TIME));

  // Following is ms.field() related, first check whether field complies with
  // LOFAR field
  if (ms.field().nrow() != 1)
    throw std::runtime_error("LOFAR MeasurementSet has multiple fields");

  if (!ms.field().tableDesc().isColumn("LOFAR_TILE_BEAM_DIR")) {
    throw std::runtime_error("LOFAR_TILE_BEAM_DIR column not found");
  }

  casacore::ScalarMeasColumn<casacore::MDirection> delay_dir_col(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));

  casacore::ArrayMeasColumn<casacore::MDirection> tile_beam_dir_col(
      ms.field(), "LOFAR_TILE_BEAM_DIR");

  // Populate struct
  ms_properties_ = {.subband_freq = band.CentreFrequency(),
                    .delay_dir = delay_dir_col(0),
                    .tile_beam_dir = *(tile_beam_dir_col(0).data())};
}

std::unique_ptr<gridded_response::GriddedResponse> LOFAR::GetGriddedResponse(
    const coords::CoordinateSystem &coordinate_system) {
  // Get and return GriddedResponse ptr
  std::unique_ptr<gridded_response::GriddedResponse> grid(
      new gridded_response::LOFARGrid(this, coordinate_system));
  // gridded_response::GriddedResponse grid(LOFARGrid(this, coordinate_system));
  return grid;
};

Station::Ptr LOFAR::ReadStation(const MeasurementSet &ms, std::size_t id,
                                const ElementResponseModel model) const {
  Station::Ptr station = ReadLofarStation(ms, id, model);
  return station;
}