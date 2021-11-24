// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "lofar.h"
#include "../griddedresponse/aartfaacgrid.h"
#include "../griddedresponse/lofargrid.h"
#include "../pointresponse/aartfaacpoint.h"
#include "../pointresponse/lofarpoint.h"
#include "../common/mathutils.h"
#include "../common/casautils.h"
#include "../msreadutils.h"
#include "../load.h"

#include <aocommon/banddata.h>
#include <cassert>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using everybeam::CorrectionMode;
using everybeam::ParseCorrectionMode;
using everybeam::Station;
using everybeam::TelescopeType;
using everybeam::ToString;
using everybeam::griddedresponse::AartfaacGrid;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::LOFARGrid;
using everybeam::pointresponse::AartfaacPoint;
using everybeam::pointresponse::LOFARPoint;
using everybeam::pointresponse::PointResponse;
using everybeam::telescope::LOFAR;

LOFAR::LOFAR(const casacore::MeasurementSet &ms, const Options &options)
    : PhasedArray(ms, options) {
  const TelescopeType telescope_type = GetTelescopeType(ms);
  if (telescope_type == TelescopeType::kAARTFAAC) {
    is_aartfaac_ = true;

    const casacore::ScalarColumn<casacore::String> antenna_type_column(
        ms.observation(), everybeam::kAartfaacAntennaTypeName);
    const std::string antenna_type = antenna_type_column(0);

    if (antenna_type != "LBA") {
      throw std::runtime_error(
          "Currently, AARTFAAC is only supported for LBA observations");
    }

    switch (options_.element_response_model) {
      case kDefault:
      case kHamaker:
        options_.element_response_model = kHamakerLba;
        break;
      case kHamakerLba:
        break;
      default:
        throw std::runtime_error(
            "Selected element response model not supported for AARTFAAC");
    }
  } else {
    if (options_.element_response_model == kDefault) {
      options_.element_response_model = kHamaker;
    }
  }

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
  if (ms.field().nrow() != 1) {
    throw std::runtime_error("LOFAR MeasurementSet has multiple fields");
  }

  if (!is_aartfaac_) {
    if (!ms.field().tableDesc().isColumn("LOFAR_TILE_BEAM_DIR")) {
      throw std::runtime_error("LOFAR_TILE_BEAM_DIR column not found");
    }
  }

  casacore::ScalarMeasColumn<casacore::MDirection> delay_dir_col(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));

  casacore::ScalarMeasColumn<casacore::MDirection> reference_dir_col(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::REFERENCE_DIR));

  CorrectionMode preapplied_correction_mode;
  casacore::MDirection preapplied_beam_dir;
  PhasedArray::CalculatePreappliedBeamOptions(ms, options_.data_column_name,
                                              preapplied_beam_dir,
                                              preapplied_correction_mode);

  size_t channel_count = band.ChannelCount();
  std::vector<double> channel_freqs(channel_count);
  for (size_t idx = 0; idx < channel_count; ++idx) {
    channel_freqs[idx] = band.ChannelFrequency(idx);
  }
  // Populate struct
  ms_properties_ = MSProperties();

  if (is_aartfaac_) {
    // Just fill with arbitrary value for AARTFAAC
    ms_properties_.tile_beam_dir = reference_dir_col(0);
  } else {
    casacore::ArrayMeasColumn<casacore::MDirection> tile_beam_dir_col(
        ms.field(), "LOFAR_TILE_BEAM_DIR");
    ms_properties_.tile_beam_dir = *(tile_beam_dir_col(0).data());
  }

  ms_properties_.subband_freq = band.ReferenceFrequency();
  ms_properties_.delay_dir = delay_dir_col(0);
  ms_properties_.reference_dir = reference_dir_col(0);
  ms_properties_.preapplied_beam_dir = preapplied_beam_dir;
  ms_properties_.preapplied_correction_mode = preapplied_correction_mode;
  ms_properties_.channel_count = channel_count;
  ms_properties_.channel_freqs = channel_freqs;
}

std::unique_ptr<GriddedResponse> LOFAR::GetGriddedResponse(
    const coords::CoordinateSystem &coordinate_system) const {
  std::unique_ptr<GriddedResponse> grid_response =
      is_aartfaac_ ? std::unique_ptr<GriddedResponse>(
                         new AartfaacGrid(this, coordinate_system))
                   : std::unique_ptr<GriddedResponse>(
                         new LOFARGrid(this, coordinate_system));
  return grid_response;
}

std::unique_ptr<PointResponse> LOFAR::GetPointResponse(double time) const {
  std::unique_ptr<PointResponse> point_response =
      is_aartfaac_
          ? std::unique_ptr<PointResponse>(new AartfaacPoint(this, time))
          : std::unique_ptr<PointResponse>(new LOFARPoint(this, time));
  return point_response;
}
