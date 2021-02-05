// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "phasedarraypoint.h"
#include "../telescope/phasedarray.h"
#include "../common/types.h"

#include "./../coords/coordutils.h"
#include "./../coords/itrfdirection.h"
#include "./../coords/itrfconverter.h"

namespace everybeam {
namespace pointresponse {

PhasedArrayPoint::PhasedArrayPoint(const telescope::Telescope* telescope_ptr,
                                   double time)
    : PointResponse(telescope_ptr, time),
      use_channel_frequency_(true),
      subband_frequency_(0.0) {
  use_differential_beam_ = telescope_->GetOptions().use_differential_beam;
}

void PhasedArrayPoint::CalculateStation(std::complex<float>* buffer, double ra,
                                        double dec, double freq,
                                        size_t station_idx, size_t field_id) {
  const telescope::PhasedArray& phasedarraytelescope =
      static_cast<const telescope::PhasedArray&>(*telescope_);

  // lock, since casacore::Direction not thread-safe
  // The lock prevents different MWAPoints to calculate the
  // the station response simultaneously
  std::unique_lock<std::mutex> lock(mtx_);
  SetITRFVectors(ra, dec);
  lock.unlock();

  double sb_freq = use_channel_frequency_ ? freq : subband_frequency_;

  if (use_differential_beam_) {
    matrix22c_t gain_matrix = phasedarraytelescope.GetStation(station_idx)
                                  ->Response(time_, freq, diff_beam_centre_,
                                             sb_freq, station0_, tile0_);
    inverse_central_gain_[0] = gain_matrix[0][0];
    inverse_central_gain_[1] = gain_matrix[0][1];
    inverse_central_gain_[2] = gain_matrix[1][0];
    inverse_central_gain_[3] = gain_matrix[1][1];
    if (!inverse_central_gain_.Invert()) {
      inverse_central_gain_ = aocommon::MC2x2F::Zero();
    }
  }

  matrix22c_t gain_matrix =
      phasedarraytelescope.GetStation(station_idx)
          ->Response(time_, freq, dir_itrf_, sb_freq, station0_, tile0_);

  if (use_differential_beam_) {
    aocommon::MC2x2F station_gains;
    station_gains[0] = gain_matrix[0][0];
    station_gains[1] = gain_matrix[0][1];
    station_gains[2] = gain_matrix[1][0];
    station_gains[3] = gain_matrix[1][1];
    aocommon::MC2x2F::ATimesB(buffer, inverse_central_gain_, station_gains);
  } else {
    buffer[0] = gain_matrix[0][0];
    buffer[1] = gain_matrix[0][1];
    buffer[2] = gain_matrix[1][0];
    buffer[3] = gain_matrix[1][1];
  }
}

void PhasedArrayPoint::SetITRFVectors(double ra, double dec) {
  coords::ITRFConverter itrf_converter(time_);
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

}  // namespace pointresponse
}  // namespace everybeam
