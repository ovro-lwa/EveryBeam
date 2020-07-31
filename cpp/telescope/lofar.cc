#include "lofar.h"
#include "../griddedresponse/lofargrid.h"
#include "../common/mathutils.h"
#include "../common/casautils.h"
#include "../lofarreadutils.h"

#include <aocommon/banddata.h>
#include <cassert>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using everybeam::Station;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::LOFARGrid;
using everybeam::telescope::LOFAR;

namespace {
bool GetPreappliedBeamDirection(casacore::MeasurementSet &ms,
                                const std::string &data_column_name,
                                bool use_differential_beam,
                                casacore::MDirection &preapplied_beam_dir) {
  casacore::ScalarMeasColumn<casacore::MDirection> referenceDirColumn(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::REFERENCE_DIR));
  preapplied_beam_dir = referenceDirColumn(0);

  // Read beam keywords of input datacolumn
  casacore::ArrayColumn<std::complex<float>> dataCol(ms, data_column_name);
  bool was_beam_applied = false;
  if (dataCol.keywordSet().isDefined("LOFAR_APPLIED_BEAM_MODE")) {
    std::string mode = dataCol.keywordSet().asString("LOFAR_APPLIED_BEAM_MODE");
    if (mode == "None")
      was_beam_applied = false;
    else {
      if (mode == "Element" || mode == "ArrayFactor")
        throw std::runtime_error(
            "This observation was corrected for the " + mode +
            " beam. WSClean can only handle a full pre-applied beam (both "
            "arrayfactor + element).");
      else if (mode == "Full") {
        was_beam_applied = true;
        casacore::String error;
        casacore::MeasureHolder mHolder;
        if (!mHolder.fromRecord(
                error, dataCol.keywordSet().asRecord("LOFAR_APPLIED_BEAM_DIR")))
          throw std::runtime_error(
              "Error while reading LOFAR_APPLIED_BEAM_DIR keyword: " + error);
        preapplied_beam_dir = mHolder.asMDirection();
      } else
        throw std::runtime_error(
            "Measurement set specifies an unknown beam correction: " + mode);
    }
  }
  if (was_beam_applied || use_differential_beam) {
    use_differential_beam = true;
  }
  return use_differential_beam;
}
}  // namespace

LOFAR::LOFAR(casacore::MeasurementSet &ms, const Options &options)
    : Telescope(ms, options) {
  stations_.resize(nstations_);
  ReadAllStations(ms, options_.element_response_model);

  // Populate MeasurementSet properties struct
  aocommon::BandData band(ms.spectralWindow());
  casacore::MSAntenna antenna(ms.antenna());
  casacore::MPosition::ScalarColumn antenna_pos_col(
      antenna, antenna.columnName(casacore::MSAntennaEnums::POSITION));
  casacore::MEpoch::ScalarColumn time_column(
      ms, ms.columnName(casacore::MSMainEnums::TIME));

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

  casacore::MDirection preapplied_beam_dir;
  options_.use_differential_beam = GetPreappliedBeamDirection(
      ms, options_.data_column_name, options_.use_differential_beam,
      preapplied_beam_dir);

  // Populate struct
  ms_properties_ = {.subband_freq = band.CentreFrequency(),
                    .delay_dir = delay_dir_col(0),
                    .tile_beam_dir = *(tile_beam_dir_col(0).data()),
                    .preapplied_beam_dir = preapplied_beam_dir};
}

std::unique_ptr<GriddedResponse> LOFAR::GetGriddedResponse(
    const coords::CoordinateSystem &coordinate_system) {
  // Get and return GriddedResponse ptr
  std::unique_ptr<GriddedResponse> grid(new LOFARGrid(this, coordinate_system));
  return grid;
};

Station::Ptr LOFAR::ReadStation(const casacore::MeasurementSet &ms,
                                std::size_t id,
                                const ElementResponseModel model) const {
  Station::Ptr station = ReadLofarStation(ms, id, model);
  return station;
}