// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mwapoint.h"
#include "../telescope/mwa.h"
#include "pointresponse/dishpoint.h"

namespace everybeam {
using mwabeam::TileBeam2016;
namespace pointresponse {

void MWAPoint::Response(BeamMode /* beam_mode */, std::complex<float>* buffer,
                        double ra, double dec, double freq,
                        size_t /* station_idx */, size_t /* field_id */) {
  const telescope::MWA& mwatelescope =
      static_cast<const telescope::MWA&>(*telescope_);

  // Only compute J2000 vectors if time was updated
  if (has_time_update_) {
    SetJ200Vectors();
    has_time_update_ = false;
  }

  if (!tile_beam_) {
    tile_beam_.reset(
        new TileBeam2016(mwatelescope.ms_properties_.delays,
                         mwatelescope.GetOptions().frequency_interpolation,
                         mwatelescope.GetOptions().coeff_path));
  }

  std::complex<double> gain[4];
  tile_beam_->ArrayResponse(ra, dec, j2000_ref_, j2000_to_hadecref_,
                            j2000_to_azelgeoref_, arr_latitude_, freq, gain);

  for (size_t i = 0; i != 4; ++i) {
    *buffer = gain[i];
    ++buffer;
  }
}

void MWAPoint::ResponseAllStations(BeamMode beam_mode,
                                   std::complex<float>* buffer, double ra,
                                   double dec, double freq, size_t) {
  Response(beam_mode, buffer, ra, dec, freq, 0.0, 0);
  // Just repeat nstations times
  for (size_t i = 1; i != telescope_->GetNrStations(); ++i) {
    std::copy_n(buffer, 4, buffer + i * 4);
  }
}

void MWAPoint::SetJ200Vectors() {
  const telescope::MWA& mwatelescope =
      static_cast<const telescope::MWA&>(*telescope_);
  // lock, since casacore::Direction not thread-safe
  // The lock prevents different MWAPoints to calculate the
  // the station response simultaneously
  std::unique_lock<std::mutex> lock(mutex_);
  casacore::MEpoch time_epoch(
      casacore::Quantity(time_ + 0.5 * update_interval_, "s"));
  casacore::MeasFrame frame(mwatelescope.ms_properties_.array_position,
                            time_epoch);

  const casacore::MDirection::Ref hadec_ref(casacore::MDirection::HADEC, frame);
  const casacore::MDirection::Ref azelgeo_ref(casacore::MDirection::AZELGEO,
                                              frame);
  j2000_ref_ = casacore::MDirection::Ref(casacore::MDirection::J2000, frame);
  j2000_to_hadecref_(j2000_ref_, hadec_ref);
  j2000_to_azelgeoref_(j2000_ref_, azelgeo_ref);
  casacore::MPosition wgs = casacore::MPosition::Convert(
      mwatelescope.ms_properties_.array_position, casacore::MPosition::WGS84)();
  arr_latitude_ = wgs.getValue().getLat();
}
}  // namespace pointresponse
}  // namespace everybeam