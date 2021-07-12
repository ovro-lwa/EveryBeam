// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "lofar.h"
#include "../griddedresponse/lofargrid.h"
#include "../pointresponse/lofarpoint.h"
#include "../common/mathutils.h"
#include "../common/casautils.h"
#include "../msreadutils.h"

#include <aocommon/banddata.h>
#include <cassert>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using everybeam::CorrectionMode;
using everybeam::ParseCorrectionMode;
using everybeam::Station;
using everybeam::ToString;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::LOFARGrid;
using everybeam::pointresponse::LOFARPoint;
using everybeam::pointresponse::PointResponse;
using everybeam::telescope::LOFAR;

namespace {
bool CalculatePreappliedBeamOptions(const casacore::MeasurementSet &ms,
                                    const std::string &data_column_name,
                                    bool force_differential_beam,
                                    casacore::MDirection &preapplied_beam_dir,
                                    CorrectionMode &correction_mode) {
  casacore::ScalarMeasColumn<casacore::MDirection> referenceDirColumn(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::REFERENCE_DIR));
  preapplied_beam_dir = referenceDirColumn(0);

  // Read beam keywords of input datacolumn
  casacore::ArrayColumn<std::complex<float>> dataCol(ms, data_column_name);
  bool was_beam_applied = false;
  if (dataCol.keywordSet().isDefined("LOFAR_APPLIED_BEAM_MODE")) {
    correction_mode = ParseCorrectionMode(
        dataCol.keywordSet().asString("LOFAR_APPLIED_BEAM_MODE"));
    switch (correction_mode) {
      case CorrectionMode::kNone:
        was_beam_applied = false;
        break;
      case CorrectionMode::kElement:
      case CorrectionMode::kArrayFactor:
      case CorrectionMode::kFull:
        was_beam_applied = true;
        casacore::String error;
        casacore::MeasureHolder mHolder;
        if (!mHolder.fromRecord(
                error, dataCol.keywordSet().asRecord("LOFAR_APPLIED_BEAM_DIR")))
          throw std::runtime_error(
              "Error while reading LOFAR_APPLIED_BEAM_DIR keyword: " + error);
        preapplied_beam_dir = mHolder.asMDirection();
        break;
    }
  } else {
    if (force_differential_beam)
      correction_mode = CorrectionMode::kFull;
    else
      correction_mode = CorrectionMode::kNone;
  }
  return was_beam_applied || force_differential_beam;
}
}  // namespace

LOFAR::LOFAR(const casacore::MeasurementSet &ms, const Options &options)
    : PhasedArray(ms, options) {
  if (options_.element_response_model == kDefault)
    options_.element_response_model = kHamaker;
  ReadAllStations(ms, stations_.begin(), options_);

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

  CorrectionMode preapplied_correction_mode;
  casacore::MDirection preapplied_beam_dir;
  options_.use_differential_beam = CalculatePreappliedBeamOptions(
      ms, options_.data_column_name, options_.use_differential_beam,
      preapplied_beam_dir, preapplied_correction_mode);

  size_t channel_count = band.ChannelCount();
  std::vector<double> channel_freqs(channel_count);
  for (size_t idx = 0; idx < channel_count; ++idx) {
    channel_freqs[idx] = band.ChannelFrequency(idx);
  }
  // Populate struct
  ms_properties_ = MSProperties();
  ms_properties_.subband_freq = band.CentreFrequency();
  ms_properties_.delay_dir = delay_dir_col(0);
  ms_properties_.tile_beam_dir = *(tile_beam_dir_col(0).data());
  ms_properties_.preapplied_beam_dir = preapplied_beam_dir;
  ms_properties_.preapplied_correction_mode = preapplied_correction_mode;
  ms_properties_.channel_count = channel_count;
  ms_properties_.channel_freqs = channel_freqs;
}

std::unique_ptr<GriddedResponse> LOFAR::GetGriddedResponse(
    const coords::CoordinateSystem &coordinate_system) const {
  // Get and return GriddedResponse ptr
  std::unique_ptr<GriddedResponse> grid(new LOFARGrid(this, coordinate_system));
  return grid;
}

std::unique_ptr<PointResponse> LOFAR::GetPointResponse(double time) const {
  std::unique_ptr<PointResponse> point_response(new LOFARPoint(this, time));
  return point_response;
}
