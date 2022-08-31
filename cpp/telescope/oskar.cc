// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oskar.h"
#include "../griddedresponse/oskargrid.h"
#include "../pointresponse/oskarpoint.h"
#include "../common/mathutils.h"
#include "../common/casautils.h"
#include "../msreadutils.h"

#include <aocommon/banddata.h>
#include <cassert>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using casacore::MeasurementSet;
using everybeam::Station;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::OSKARGrid;
using everybeam::pointresponse::OSKARPoint;
using everybeam::pointresponse::PointResponse;
using everybeam::telescope::OSKAR;

OSKAR::OSKAR(const MeasurementSet& ms, const Options& options)
    : PhasedArray(ms, options) {
  if (options_.element_response_model == ElementResponseModel::kDefault) {
    options_.element_response_model = ElementResponseModel::kOSKARSphericalWave;
  }
  ReadAllStations(ms, stations_.begin(), options_);

  aocommon::BandData band(ms.spectralWindow());
  casacore::ScalarMeasColumn<casacore::MDirection> delay_dir_col(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));

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
  ms_properties_.subband_freq = band.ReferenceFrequency();
  ms_properties_.delay_dir = delay_dir_col(0);
  // tile_beam_dir has dummy values for OSKAR
  ms_properties_.tile_beam_dir = delay_dir_col(0);
  ms_properties_.preapplied_beam_dir = preapplied_beam_dir;
  ms_properties_.preapplied_correction_mode = preapplied_correction_mode;
  ms_properties_.channel_count = channel_count;
  ms_properties_.channel_freqs = channel_freqs;
}

std::unique_ptr<GriddedResponse> OSKAR::GetGriddedResponse(
    const coords::CoordinateSystem& coordinate_system) const {
  return std::make_unique<OSKARGrid>(this, coordinate_system);
}

std::unique_ptr<PointResponse> OSKAR::GetPointResponse(double time) const {
  return std::make_unique<OSKARPoint>(this, time);
}