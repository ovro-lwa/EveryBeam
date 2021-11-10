// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "phasedarraygrid.h"
#include "../telescope/phasedarray.h"

#include <aocommon/imagecoordinates.h>
#include <aocommon/threadpool.h>
#include <cmath>
#include <iostream>

using everybeam::griddedresponse::PhasedArrayGrid;

PhasedArrayGrid::PhasedArrayGrid(
    const telescope::Telescope* telescope_ptr,
    const coords::CoordinateSystem& coordinate_system)
    : GriddedResponse(telescope_ptr, coordinate_system),
      beam_normalisation_mode_(
          telescope_->GetOptions().beam_normalisation_mode) {
  // Extract PhasedArrayPoint specific options from ms_properties_ and
  // telescope::Options
  const telescope::PhasedArray& phasedarray =
      dynamic_cast<const telescope::PhasedArray&>(*telescope_ptr);
  delay_dir_ = phasedarray.GetMSProperties().delay_dir;
  tile_beam_dir_ = phasedarray.GetMSProperties().tile_beam_dir;
  preapplied_beam_dir_ = phasedarray.GetMSProperties().preapplied_beam_dir;
  preapplied_correction_mode_ =
      phasedarray.GetMSProperties().preapplied_correction_mode;
  subband_frequency_ = phasedarray.GetMSProperties().subband_freq;
  use_channel_frequency_ = phasedarray.GetOptions().use_channel_frequency;

  // Compute and set number of threads
  const size_t ncpus = aocommon::ThreadPool::NCPUs();
  const size_t nthreads = std::min(ncpus, telescope_->GetNrStations());
  threads_.resize(nthreads);
}

void PhasedArrayGrid::Response(BeamMode beam_mode, std::complex<float>* buffer,
                               double time, double frequency,
                               size_t station_idx,
                               [[maybe_unused]] size_t field_id) {
  aocommon::Lane<Job> lane(threads_.size());
  lane_ = &lane;

  SetITRFVectors(time);

  inverse_central_gain_.resize(1);
  bool apply_normalisation = CalculateBeamNormalisation(
      beam_mode, time, frequency, station_idx, inverse_central_gain_[0]);

  // Prepare threads
  for (auto& thread : threads_) {
    thread = std::thread(&PhasedArrayGrid::CalcThread, this, beam_mode,
                         apply_normalisation, buffer, time, frequency);
  }

  for (size_t y = 0; y != height_; ++y) {
    lane.emplace(Job(y, station_idx, 0));
  }

  lane.write_end();
  for (auto& thread : threads_) thread.join();
}

void PhasedArrayGrid::ResponseAllStations(BeamMode beam_mode,
                                          std::complex<float>* buffer,
                                          double time, double frequency,
                                          size_t) {
  const telescope::PhasedArray& phasedarraytelescope =
      static_cast<const telescope::PhasedArray&>(*telescope_);
  aocommon::Lane<Job> lane(threads_.size());
  lane_ = &lane;

  SetITRFVectors(time);

  bool apply_normalisation;
  inverse_central_gain_.resize(phasedarraytelescope.GetNrStations());
  for (size_t i = 0; i != phasedarraytelescope.GetNrStations(); ++i) {
    apply_normalisation = CalculateBeamNormalisation(
        beam_mode, time, frequency, i, inverse_central_gain_[i]);
  }

  // Prepare threads
  for (auto& thread : threads_) {
    thread = std::thread(&PhasedArrayGrid::CalcThread, this, beam_mode,
                         apply_normalisation, buffer, time, frequency);
  }

  for (size_t y = 0; y != height_; ++y) {
    for (size_t antenna_idx = 0;
         antenna_idx != phasedarraytelescope.GetNrStations(); ++antenna_idx) {
      lane.write(Job(y, antenna_idx, antenna_idx));
    }
  }

  lane.write_end();
  for (auto& thread : threads_) thread.join();
}

void PhasedArrayGrid::SetITRFVectors(double time) {
  coords::ITRFConverter itrf_converter(time);
  coords::SetITRFVector(itrf_converter.ToDirection(delay_dir_), station0_);
  coords::SetITRFVector(itrf_converter.ToDirection(tile_beam_dir_), tile0_);

  const casacore::Unit rad_unit("rad");

  casacore::MDirection l_dir(
      casacore::MVDirection(casacore::Quantity(ra_ + M_PI / 2, rad_unit),
                            casacore::Quantity(0, rad_unit)),
      casacore::MDirection::J2000);
  coords::SetITRFVector(itrf_converter.ToDirection(l_dir), l_vector_itrf_);

  casacore::MDirection m_dir(
      casacore::MVDirection(casacore::Quantity(ra_, rad_unit),
                            casacore::Quantity(dec_ + M_PI / 2, rad_unit)),
      casacore::MDirection::J2000);
  coords::SetITRFVector(itrf_converter.ToDirection(m_dir), m_vector_itrf_);

  casacore::MDirection n_dir(
      casacore::MVDirection(casacore::Quantity(ra_, rad_unit),
                            casacore::Quantity(dec_, rad_unit)),
      casacore::MDirection::J2000);
  coords::SetITRFVector(itrf_converter.ToDirection(n_dir), n_vector_itrf_);

  coords::SetITRFVector(itrf_converter.ToDirection(preapplied_beam_dir_),
                        diff_beam_centre_);
}

bool PhasedArrayGrid::CalculateBeamNormalisation(
    BeamMode beam_mode, double time, double frequency, size_t station_idx,
    aocommon::MC2x2F& inverse_gain) const {
  const telescope::PhasedArray& phasedarraytelescope =
      static_cast<const telescope::PhasedArray&>(*telescope_);
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
    // in all other cases the response for the reference direction with
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

void PhasedArrayGrid::CalcThread(BeamMode beam_mode, bool apply_normalisation,
                                 std::complex<float>* buffer, double time,
                                 double frequency) {
  const telescope::PhasedArray& phasedarraytelescope =
      static_cast<const telescope::PhasedArray&>(*telescope_);
  const size_t values_per_ant = width_ * height_ * 4;
  const double sb_freq =
      use_channel_frequency_ ? frequency : subband_frequency_;

  Job job;
  while (lane_->read(job)) {
    for (size_t x = 0; x != width_; ++x) {
      double l, m, n;
      aocommon::ImageCoordinates::XYToLM(x, job.y, dl_, dm_, width_, height_, l,
                                         m);
      l += phase_centre_dl_;
      m += phase_centre_dm_;
      n = sqrt(1.0 - l * l - m * m);

      const vector3r_t itrf_direction =
          l * l_vector_itrf_ + m * m_vector_itrf_ + n * n_vector_itrf_;

      std::complex<float>* base_buffer = buffer + (x + job.y * width_) * 4;

      std::complex<float>* ant_buffer_ptr =
          base_buffer + job.buffer_offset * values_per_ant;

      const aocommon::MC2x2F gain_matrix = aocommon::MC2x2F(
          phasedarraytelescope.GetStation(job.antenna_idx)
              ->Response(beam_mode, time, frequency, itrf_direction, sb_freq,
                         station0_, tile0_)
              .Data());

      if (apply_normalisation) {
        aocommon::MC2x2F::ATimesB(ant_buffer_ptr,
                                  inverse_central_gain_[job.buffer_offset],
                                  gain_matrix);
      } else {
        gain_matrix.AssignTo(ant_buffer_ptr);
      }
    }
  }
}
