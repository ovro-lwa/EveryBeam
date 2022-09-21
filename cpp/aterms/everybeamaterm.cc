// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "everybeamaterm.h"

#include "../load.h"
#include "../options.h"
#include "../elementresponse.h"
#include "../griddedresponse/griddedresponse.h"

namespace everybeam {
namespace aterms {

EveryBeamATerm::EveryBeamATerm(
    const casacore::MeasurementSet& ms,
    const aocommon::CoordinateSystem& coordinate_system, const Options& options)
    : telescope_(Load(ms, options)),
      coordinate_system_(coordinate_system),
      beam_mode_(options.beam_mode) {}

bool EveryBeamATerm::CalculateBeam(std::complex<float>* buffer, double time,
                                   double frequency, size_t field_id) {
  if (!telescope_->GetIsTimeRelevant()) {
    if (field_id == cached_field_id && cached_frequency_ == frequency) {
      // Exit calculation
      return false;
    } else {
      // Update cached values
      cached_field_id = field_id;
      cached_frequency_ = frequency;
    }
  }

  // Get the gridded response
  std::unique_ptr<griddedresponse::GriddedResponse> grid_response =
      telescope_->GetGriddedResponse(coordinate_system_);

  grid_response->ResponseAllStations(beam_mode_, buffer, time, frequency,
                                     field_id);

  SaveATermsIfNecessary(buffer, telescope_->GetNrStations(),
                        coordinate_system_.width, coordinate_system_.height);
  return true;
}

}  // namespace aterms
}  // namespace everybeam
