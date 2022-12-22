// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "phasedarraypoint.h"
#include "../telescope/phasedarray.h"
#include "../common/types.h"

#include "./../coords/itrfdirection.h"
#include "./../coords/itrfconverter.h"

#include <limits>
namespace everybeam {

using telescope::PhasedArray;

namespace pointresponse {

PhasedArrayPoint::PhasedArrayPoint(const telescope::Telescope* telescope_ptr,
                                   double time)
    : PointResponse(telescope_ptr, time),
      PhasedArrayResponse(static_cast<const PhasedArray*>(telescope_ptr)),
      ra_(std::numeric_limits<double>::min()),
      dec_(std::numeric_limits<double>::min()),
      has_partial_itrf_update_(false),
      is_local_(false),
      rotate_(true) {}

void PhasedArrayPoint::Response(BeamMode beam_mode, std::complex<float>* buffer,
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

  const PhasedArray& phased_array =
      static_cast<const PhasedArray&>(*telescope_);

  aocommon::MC2x2F inverse_central_gain;
  const bool apply_normalisation = CalculateBeamNormalisation(
      beam_mode, time_, freq, station_idx, inverse_central_gain);
  const double sb_freq = use_channel_frequency_ ? freq : subband_frequency_;

  const aocommon::MC2x2F gain_matrix =
      aocommon::MC2x2F(phased_array.GetStation(station_idx)
                           .Response(beam_mode, time_, freq, itrf_direction_,
                                     sb_freq, station0_, tile0_)
                           .Data());

  if (apply_normalisation) {
    aocommon::MC2x2F::ATimesB(buffer, inverse_central_gain, gain_matrix);
  } else {
    gain_matrix.AssignTo(buffer);
  }
}

aocommon::MC2x2 PhasedArrayPoint::Response(BeamMode beam_mode,
                                           size_t station_idx, double freq,
                                           const vector3r_t& direction,
                                           std::mutex* mutex) {
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

  aocommon::MC2x2F inverse_central_gain;
  const bool apply_normalisation = CalculateBeamNormalisation(
      beam_mode, time_, freq, station_idx, inverse_central_gain);

  aocommon::MC2x2 gain_matrix = UnnormalisedResponse(
      beam_mode, station_idx, freq, direction, station0_, tile0_);

  // Conversion to MC2x2 (double) for inverse_central_gain needed
  return apply_normalisation
             ? aocommon::MC2x2(inverse_central_gain.Data()) * gain_matrix
             : gain_matrix;
}

aocommon::MC2x2 PhasedArrayPoint::UnnormalisedResponse(
    BeamMode beam_mode, size_t station_idx, double freq,
    const vector3r_t& direction, const vector3r_t& station0,
    const vector3r_t& tile0) const {
  const PhasedArray& phased_array =
      static_cast<const PhasedArray&>(*telescope_);
  const double sb_freq = use_channel_frequency_ ? freq : subband_frequency_;

  return phased_array.GetStation(station_idx)
      .Response(beam_mode, time_, freq, direction, sb_freq, station0, tile0,
                rotate_);
}

aocommon::MC2x2 PhasedArrayPoint::ElementResponse(size_t station_idx,
                                                  double freq,
                                                  const vector3r_t& direction,
                                                  size_t element_idx) const {
  const PhasedArray& phased_array =
      static_cast<const PhasedArray&>(*telescope_);
  return phased_array.GetStation(station_idx)
      .ComputeElementResponse(time_, freq, direction, element_idx, is_local_,
                              rotate_);
}

void PhasedArrayPoint::UpdateITRFVectors(double ra, double dec) {
  ra_ = ra;
  dec_ = dec;
  // lock, since casacore::Direction is not thread-safe
  // The lock prevents different PhasedArrayPoints to calculate the
  // the station response simultaneously
  std::unique_lock<std::mutex> lock(mutex_);
  const coords::ItrfConverter itrf_converter(time_ + 0.5 * update_interval_);
  station0_ = itrf_converter.ToItrf(delay_dir_);
  tile0_ = itrf_converter.ToItrf(tile_beam_dir_);
  // Only the n vector is relevant for a single point. l and m are not.
  itrf_direction_ = itrf_converter.RaDecToItrf(ra, dec);
  diff_beam_centre_ = itrf_converter.ToItrf(preapplied_beam_dir_);
}

void PhasedArrayPoint::UpdateITRFVectors(std::mutex& mutex) {
  std::unique_lock<std::mutex> lock(mutex);
  const coords::ItrfConverter itrf_converter(time_);
  station0_ = itrf_converter.ToItrf(delay_dir_);
  tile0_ = itrf_converter.ToItrf(tile_beam_dir_);
}

}  // namespace pointresponse
}  // namespace everybeam
