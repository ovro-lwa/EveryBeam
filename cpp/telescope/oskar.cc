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

OSKAR::OSKAR(const MeasurementSet &ms, const Options &options)
    : PhasedArray(ms, options) {
  if (options_.element_response_model == kDefault)
    options_.element_response_model = kOSKARSphericalWave;
  ReadAllStations(ms, options_.element_response_model);

  aocommon::BandData band(ms.spectralWindow());
  casacore::ScalarMeasColumn<casacore::MDirection> delay_dir_col(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));

  size_t channel_count = band.ChannelCount();
  std::vector<double> channel_freqs(channel_count);
  for (size_t idx = 0; idx < channel_count; ++idx) {
    channel_freqs[idx] = band.ChannelFrequency(idx);
  }

  // Populate struct
  ms_properties_ = MSProperties();
  ms_properties_.subband_freq = band.CentreFrequency();
  ms_properties_.delay_dir = delay_dir_col(0);
  // tile_beam_dir and preapplied_beam_dir
  // have dummy values for OSKAR
  ms_properties_.tile_beam_dir = delay_dir_col(0);
  ms_properties_.preapplied_beam_dir = delay_dir_col(0);
  ms_properties_.channel_count = channel_count;
  ms_properties_.channel_freqs = channel_freqs;
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
