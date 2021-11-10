// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "phasedarraypoint.h"
#include "../telescope/phasedarray.h"
#include "../common/types.h"

#include "./../coords/coordutils.h"
#include "./../coords/itrfdirection.h"
#include "./../coords/itrfconverter.h"

#include <limits>
namespace everybeam {
namespace pointresponse {

PhasedArrayPoint::PhasedArrayPoint(const telescope::Telescope *telescope_ptr,
                                   double time)
    : PointResponse(telescope_ptr, time),
      beam_normalisation_mode_(
          telescope_->GetOptions().beam_normalisation_mode),
      ra_(std::numeric_limits<double>::min()),
      dec_(std::numeric_limits<double>::min()),
      has_partial_itrf_update_(false),
      is_local_(false),
      rotate_(true) {
  // Extract phased array specific options from ms_properties_ and
  // telescope::Options
  // TODO: code is largely a duplicate of PhasedArrayGrid
  // constructor
  const telescope::PhasedArray &phasedarray =
      dynamic_cast<const telescope::PhasedArray &>(*telescope_ptr);
  delay_dir_ = phasedarray.GetMSProperties().delay_dir;
  tile_beam_dir_ = phasedarray.GetMSProperties().tile_beam_dir;
  preapplied_beam_dir_ = phasedarray.GetMSProperties().preapplied_beam_dir;
  preapplied_correction_mode_ =
      phasedarray.GetMSProperties().preapplied_correction_mode;
  subband_frequency_ = phasedarray.GetMSProperties().subband_freq;
  use_channel_frequency_ = phasedarray.GetOptions().use_channel_frequency;
}

void PhasedArrayPoint::Response(BeamMode beam_mode, std::complex<float> *buffer,
                                double ra, double dec, double freq,
                                size_t station_idx,
                                [[maybe_unused]] size_t field_id) {
  // Only compute ITRF directions if values differ from cached values
  if (has_time_update_ || has_partial_itrf_update_ ||
      std::abs(ra - ra_) > 1e-10 || std::abs(dec - dec_) > 1e-10) {
    UpdateITRFVectors(ra, dec);
    has_time_update_ = false;
    has_partial_itrf_update_ = false;
  }

  const telescope::PhasedArray &phasedarraytelescope =
      static_cast<const telescope::PhasedArray &>(*telescope_);

  aocommon::MC2x2F inverse_central_gain;
  bool apply_normalisation = CalculateBeamNormalisation(
      beam_mode, time_, freq, station_idx, inverse_central_gain);
  const double sb_freq = use_channel_frequency_ ? freq : subband_frequency_;

  const aocommon::MC2x2F gain_matrix =
      aocommon::MC2x2F(phasedarraytelescope.GetStation(station_idx)
                           ->Response(beam_mode, time_, freq, dir_itrf_,
                                      sb_freq, station0_, tile0_)
                           .Data());

  if (apply_normalisation) {
    aocommon::MC2x2F::ATimesB(buffer, inverse_central_gain, gain_matrix);
  } else {
    gain_matrix.AssignTo(buffer);
  }
}

aocommon::MC2x2 PhasedArrayPoint::FullResponse(size_t station_idx, double freq,
                                               const vector3r_t &direction,
                                               std::mutex *mutex) {
  if (has_time_update_) {
    if (mutex != nullptr) {
      // Caller takes over responsibility to be thread-safe
      UpdateITRFVectors(*mutex);
    } else {
      // Callee assumes that caller is thread-safe
      UpdateITRFVectors(mutex_);
    }
    has_time_update_ = false;
    has_partial_itrf_update_ = true;
  }
  return FullResponse(station_idx, freq, direction, station0_, tile0_);
}

aocommon::MC2x2 PhasedArrayPoint::FullResponse(size_t station_idx, double freq,
                                               const vector3r_t &direction,
                                               const vector3r_t &station0,
                                               const vector3r_t &tile0) {
  const telescope::PhasedArray &phasedarraytelescope =
      static_cast<const telescope::PhasedArray &>(*telescope_);
  const double sb_freq = use_channel_frequency_ ? freq : subband_frequency_;
  return phasedarraytelescope.GetStation(station_idx)
      ->Response(time_, freq, direction, sb_freq, station0, tile0, rotate_);
}

aocommon::MC2x2Diag PhasedArrayPoint::ArrayFactor(size_t station_idx,
                                                  double freq,
                                                  const vector3r_t &direction,
                                                  std::mutex *mutex) {
  if (has_time_update_) {
    if (mutex != nullptr) {
      // Caller takes over responsibility to be thread-safe
      UpdateITRFVectors(*mutex);
    } else {
      // Calllee assumes that caller is thread-safe
      UpdateITRFVectors(mutex_);
    }
    has_time_update_ = false;
    has_partial_itrf_update_ = true;
  }
  return ArrayFactor(station_idx, freq, direction, station0_, tile0_);
}

aocommon::MC2x2Diag PhasedArrayPoint::ArrayFactor(size_t station_idx,
                                                  double freq,
                                                  const vector3r_t &direction,
                                                  const vector3r_t &station0,
                                                  const vector3r_t &tile0) {
  const telescope::PhasedArray &phasedarraytelescope =
      static_cast<const telescope::PhasedArray &>(*telescope_);

  const double sb_freq = use_channel_frequency_ ? freq : subband_frequency_;

  return phasedarraytelescope.GetStation(station_idx)
      ->ArrayFactor(time_, freq, direction, sb_freq, station0, tile0);
}

aocommon::MC2x2 PhasedArrayPoint::ElementResponse(
    size_t station_idx, double freq, const vector3r_t &direction) const {
  const telescope::PhasedArray &phasedarraytelescope =
      static_cast<const telescope::PhasedArray &>(*telescope_);
  return phasedarraytelescope.GetStation(station_idx)
      ->ComputeElementResponse(time_, freq, direction, is_local_, rotate_);
}

aocommon::MC2x2 PhasedArrayPoint::ElementResponse(size_t station_idx,
                                                  double freq,
                                                  const vector3r_t &direction,
                                                  size_t element_idx) const {
  const telescope::PhasedArray &phasedarraytelescope =
      static_cast<const telescope::PhasedArray &>(*telescope_);
  return phasedarraytelescope.GetStation(station_idx)
      ->ComputeElementResponse(time_, freq, direction, element_idx, is_local_,
                               rotate_);
}

void PhasedArrayPoint::UpdateITRFVectors(double ra, double dec) {
  ra_ = ra;
  dec_ = dec;
  // lock, since casacore::Direction is not thread-safe
  // The lock prevents different PhasedArrayPoints to calculate the
  // the station response simultaneously
  std::unique_lock<std::mutex> lock(mutex_);
  coords::ITRFConverter itrf_converter(time_ + 0.5 * update_interval_);
  coords::SetITRFVector(itrf_converter.ToDirection(delay_dir_), station0_);
  coords::SetITRFVector(itrf_converter.ToDirection(tile_beam_dir_), tile0_);

  const casacore::Unit rad_unit("rad");

  // Only n_dir relevant for a single point
  casacore::MDirection n_dir(
      casacore::MVDirection(casacore::Quantity(ra, rad_unit),
                            casacore::Quantity(dec, rad_unit)),
      casacore::MDirection::J2000);
  coords::SetITRFVector(itrf_converter.ToDirection(n_dir), dir_itrf_);

  coords::SetITRFVector(itrf_converter.ToDirection(preapplied_beam_dir_),
                        diff_beam_centre_);
}

void PhasedArrayPoint::UpdateITRFVectors(std::mutex &mutex) {
  std::unique_lock<std::mutex> lock(mutex);
  coords::ITRFConverter itrf_converter(time_);
  coords::SetITRFVector(itrf_converter.ToDirection(delay_dir_), station0_);
  coords::SetITRFVector(itrf_converter.ToDirection(tile_beam_dir_), tile0_);
}

bool PhasedArrayPoint::CalculateBeamNormalisation(
    BeamMode beam_mode, double time, double frequency, size_t station_idx,
    aocommon::MC2x2F &inverse_gain) const {
  const telescope::PhasedArray &phasedarraytelescope =
      static_cast<const telescope::PhasedArray &>(*telescope_);
  if (beam_normalisation_mode_ == BeamNormalisationMode::kNone) {
    return false;
  }

  const double sb_freq =
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
        phasedarraytelescope.GetStation(station_idx)
            ->Response(preapplied_correction_mode_, time, frequency,
                       diff_beam_centre_, sb_freq, station0_, tile0_)
            .Data());
  } else {
    // in all other cases the response in for the reference direction with
    // beam_mode is needed
    inverse_gain = aocommon::MC2x2F(
        phasedarraytelescope.GetStation(station_idx)
            ->Response(beam_mode, time, frequency, diff_beam_centre_, sb_freq,
                       station0_, tile0_)
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
      const float amplitude_inv = 1.0 / std::sqrt(0.5 * Norm(inverse_gain));
      inverse_gain[0] = std::isfinite(amplitude_inv) ? amplitude_inv : 0.0;
      inverse_gain[1] = 0.0;
      inverse_gain[2] = 0.0;
      inverse_gain[3] = std::isfinite(amplitude_inv) ? amplitude_inv : 0.0;
      break;
    }
    default:
      throw std::runtime_error("Invalid beam normalisation mode here");
  }

  return true;
}

}  // namespace pointresponse
}  // namespace everybeam
