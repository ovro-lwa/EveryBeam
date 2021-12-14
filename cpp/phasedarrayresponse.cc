// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "telescope/phasedarray.h"
#include "phasedarrayresponse.h"

#include <cassert>

namespace everybeam {

PhasedArrayResponse::PhasedArrayResponse(
    const telescope::PhasedArray* phased_array)
    : delay_dir_(phased_array->GetMSProperties().delay_dir),
      tile_beam_dir_(phased_array->GetMSProperties().tile_beam_dir),
      preapplied_beam_dir_(phased_array->GetMSProperties().preapplied_beam_dir),
      preapplied_correction_mode_(
          phased_array->GetMSProperties().preapplied_correction_mode),
      beam_normalisation_mode_(
          phased_array->GetOptions().beam_normalisation_mode),
      use_channel_frequency_(phased_array->GetOptions().use_channel_frequency),
      subband_frequency_(phased_array->GetMSProperties().subband_freq),
      phased_array_(phased_array) {
  assert(phased_array_);
}

bool PhasedArrayResponse::CalculateBeamNormalisation(
    BeamMode beam_mode, double time, double frequency, size_t station_index,
    aocommon::MC2x2F& inverse_gain) const {
  const telescope::PhasedArray& phased_array =
      static_cast<const telescope::PhasedArray&>(*phased_array_);
  if (beam_normalisation_mode_ == BeamNormalisationMode::kNone) {
    return false;
  }

  const double subband_frequency =
      use_channel_frequency_ ? frequency : subband_frequency_;

  // if the normalisation mode is kPreApplied, but no beam correction was pre
  // applied then there is nothing to do
  if (beam_normalisation_mode_ == BeamNormalisationMode::kPreApplied &&
      preapplied_correction_mode_ == BeamMode::kNone) {
    return false;
  }

  // If the normalisation mode is kPreApplied, or kPreAppliedOrFull and the
  // fallback to Full is not needed then the response for the diff_beam_centre_
  // with preapplied_correction_mode_ needs to be computed
  if (beam_normalisation_mode_ == BeamNormalisationMode::kPreApplied ||
      (beam_normalisation_mode_ == BeamNormalisationMode::kPreAppliedOrFull &&
       preapplied_correction_mode_ != BeamMode::kNone)) {
    inverse_gain = aocommon::MC2x2F(
        phased_array.GetStation(station_index)
            ->Response(preapplied_correction_mode_, time, frequency,
                       diff_beam_centre_, subband_frequency, station0_, tile0_)
            .Data());
  } else {
    // in all other cases the response for the reference direction with
    // beam_mode is needed
    inverse_gain = aocommon::MC2x2F(
        phased_array.GetStation(station_index)
            ->Response(beam_mode, time, frequency, diff_beam_centre_,
                       subband_frequency, station0_, tile0_)
            .Data());
  }

  switch (beam_normalisation_mode_) {
    case BeamNormalisationMode::kFull:
    case BeamNormalisationMode::kPreApplied:
    case BeamNormalisationMode::kPreAppliedOrFull:
      if (!inverse_gain.Invert()) {
        inverse_gain = aocommon::MC2x2F::Zero();
      }
      break;
    case BeamNormalisationMode::kAmplitude: {
      const float norm_inverse_gain = Norm(inverse_gain);
      const float amplitude_inv =
          (norm_inverse_gain == 0.0) ? 0.0
                                     : 1.0 / std::sqrt(0.5 * norm_inverse_gain);
      inverse_gain[0] = amplitude_inv;
      inverse_gain[1] = 0.0;
      inverse_gain[2] = 0.0;
      inverse_gain[3] = amplitude_inv;
      break;
    }
    case BeamNormalisationMode::kNone:
      throw std::runtime_error("Invalid beam normalisation mode here");
  }
  return true;
}

}  // namespace everybeam