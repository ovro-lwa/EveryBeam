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
      use_channel_frequency_(telescope_->GetOptions().use_channel_frequency),
      use_differential_beam_(telescope_->GetOptions().use_differential_beam),
      preapplied_correction_mode_(CorrectionMode::kFull),
      subband_frequency_(0.0),
      ra_(std::numeric_limits<double>::min()),
      dec_(std::numeric_limits<double>::min()),
      has_partial_itrf_update_(false),
      is_local_(false),
      rotate_(true) {}

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

  const double sb_freq = use_channel_frequency_ ? freq : subband_frequency_;

  aocommon::MC2x2F inverse_central_gain;
  if (use_differential_beam_) {
    inverse_central_gain = aocommon::MC2x2F(
        phasedarraytelescope.GetStation(station_idx)
            ->Response(preapplied_correction_mode_, time_, freq,
                       diff_beam_centre_, sb_freq, station0_, tile0_)
            .Data());
    if (!inverse_central_gain.Invert()) {
      inverse_central_gain = aocommon::MC2x2F::Zero();
    }
  }

  // TODO: alternative could be a call to FullResponse()
  const aocommon::MC2x2F gain_matrix = aocommon::MC2x2F(
      phasedarraytelescope.GetStation(station_idx)
          ->Response(time_, freq, dir_itrf_, sb_freq, station0_, tile0_)
          .Data());

  if (use_differential_beam_) {
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
}  // namespace pointresponse
}  // namespace everybeam
