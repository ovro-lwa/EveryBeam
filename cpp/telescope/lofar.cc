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

namespace {
bool GetPreappliedBeamDirection(casacore::MeasurementSet &ms,
                                const std::string &dataColumnName,
                                bool useDifferentialBeam,
                                casacore::MDirection &preappliedBeamDir) {
  casacore::ScalarMeasColumn<casacore::MDirection> referenceDirColumn(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::REFERENCE_DIR));
  preappliedBeamDir = referenceDirColumn(0);

  // Read beam keywords of input datacolumn
  casacore::ArrayColumn<std::complex<float>> dataCol(ms, dataColumnName);
  bool wasBeamApplied = false;
  if (dataCol.keywordSet().isDefined("LOFAR_APPLIED_BEAM_MODE")) {
    std::string mode = dataCol.keywordSet().asString("LOFAR_APPLIED_BEAM_MODE");
    if (mode == "None")
      wasBeamApplied = false;
    else {
      if (mode == "Element" || mode == "ArrayFactor")
        throw std::runtime_error(
            "This observation was corrected for the " + mode +
            " beam. WSClean can only handle a full pre-applied beam (both "
            "arrayfactor + element).");
      else if (mode == "Full") {
        wasBeamApplied = true;
        casacore::String error;
        casacore::MeasureHolder mHolder;
        if (!mHolder.fromRecord(
                error, dataCol.keywordSet().asRecord("LOFAR_APPLIED_BEAM_DIR")))
          throw std::runtime_error(
              "Error while reading LOFAR_APPLIED_BEAM_DIR keyword: " + error);
        preappliedBeamDir = mHolder.asMDirection();
      } else
        throw std::runtime_error(
            "Measurement set specifies an unknown beam correction: " + mode);
    }
  }
  if (wasBeamApplied || useDifferentialBeam) {
    useDifferentialBeam = true;
  }
  return useDifferentialBeam;
}
}  // namespace

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