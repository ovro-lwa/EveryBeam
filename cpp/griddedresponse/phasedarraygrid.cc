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
      use_channel_frequency_(true),
      subband_frequency_(0.0) {
  // Set private members
  use_differential_beam_ = telescope_->GetOptions().use_differential_beam;
  size_t ncpus = aocommon::ThreadPool::NCPUs();
  nthreads_ = std::min(ncpus, telescope_->GetNrStations());
  threads_.resize(nthreads_);
}

void PhasedArrayGrid::CalculateStation(std::complex<float>* buffer, double time,
                                       double frequency, size_t station_idx,
                                       size_t field_id) {
  const telescope::PhasedArray& phasedarraytelescope =
      static_cast<const telescope::PhasedArray&>(*telescope_);
  aocommon::Lane<Job> lane(nthreads_);
  lane_ = &lane;

  SetITRFVectors(time);
  if (use_differential_beam_) {
    double sb_freq = use_channel_frequency_ ? frequency : subband_frequency_;
    inverse_central_gain_.resize(1);
    inverse_central_gain_[0] = aocommon::MC2x2F(
        phasedarraytelescope.GetStation(station_idx)
            ->Response(preapplied_correction_mode_, time, frequency,
                       diff_beam_centre_, sb_freq, station0_, tile0_)
            .Data());
    if (!inverse_central_gain_[0].Invert()) {
      inverse_central_gain_[0] = aocommon::MC2x2F::Zero();
    }
  }

  // Prepare threads
  for (size_t i = 0; i != nthreads_; ++i) {
    threads_[i] = std::thread(&PhasedArrayGrid::CalcThread, this, buffer, time,
                              frequency);
  }

  for (size_t y = 0; y != height_; ++y) {
    lane.emplace(Job(y, station_idx, 0));
  }

  lane.write_end();
  for (size_t i = 0; i != nthreads_; ++i) threads_[i].join();
}

void PhasedArrayGrid::CalculateAllStations(std::complex<float>* buffer,
                                           double time, double frequency,
                                           size_t) {
  const telescope::PhasedArray& phasedarraytelescope =
      static_cast<const telescope::PhasedArray&>(*telescope_);
  aocommon::Lane<Job> lane(nthreads_);
  lane_ = &lane;

  SetITRFVectors(time);

  if (use_differential_beam_) {
    double sb_freq = use_channel_frequency_ ? frequency : subband_frequency_;
    inverse_central_gain_.resize(phasedarraytelescope.GetNrStations());
    for (size_t i = 0; i != phasedarraytelescope.GetNrStations(); ++i) {
      inverse_central_gain_[i] =
          aocommon::MC2x2F(phasedarraytelescope.GetStation(i)
                               ->Response(time, frequency, diff_beam_centre_,
                                          sb_freq, station0_, tile0_)
                               .Data());
      if (!inverse_central_gain_[i].Invert()) {
        inverse_central_gain_[i] = aocommon::MC2x2F::Zero();
      }
    }
  }

  // Prepare threads
  for (size_t i = 0; i != nthreads_; ++i) {
    threads_[i] = std::thread(&PhasedArrayGrid::CalcThread, this, buffer, time,
                              frequency);
  }

  for (size_t y = 0; y != height_; ++y) {
    for (size_t antenna_idx = 0;
         antenna_idx != phasedarraytelescope.GetNrStations(); ++antenna_idx) {
      lane.write(Job(y, antenna_idx, antenna_idx));
    }
  }

  lane.write_end();
  for (size_t i = 0; i != nthreads_; ++i) threads_[i].join();
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

void PhasedArrayGrid::CalcThread(std::complex<float>* buffer, double time,
                                 double frequency) {
  const telescope::PhasedArray& phasedarraytelescope =
      static_cast<const telescope::PhasedArray&>(*telescope_);
  const size_t values_per_ant = width_ * height_ * 4;
  double sb_freq = use_channel_frequency_ ? frequency : subband_frequency_;

  Job job;
  while (lane_->read(job)) {
    for (size_t x = 0; x != width_; ++x) {
      double l, m, n;
      aocommon::ImageCoordinates::XYToLM(x, job.y, dl_, dm_, width_, height_, l,
                                         m);
      l += phase_centre_dl_;
      m += phase_centre_dm_;
      n = sqrt(1.0 - l * l - m * m);

      vector3r_t itrf_direction;

      itrf_direction[0] =
          l * l_vector_itrf_[0] + m * m_vector_itrf_[0] + n * n_vector_itrf_[0];
      itrf_direction[1] =
          l * l_vector_itrf_[1] + m * m_vector_itrf_[1] + n * n_vector_itrf_[1];
      itrf_direction[2] =
          l * l_vector_itrf_[2] + m * m_vector_itrf_[2] + n * n_vector_itrf_[2];

      std::complex<float>* base_buffer = buffer + (x + job.y * width_) * 4;

      std::complex<float>* ant_buffer_ptr =
          base_buffer + job.buffer_offset * values_per_ant;

      const aocommon::MC2x2F gain_matrix =
          aocommon::MC2x2F(phasedarraytelescope.GetStation(job.antenna_idx)
                               ->Response(time, frequency, itrf_direction,
                                          sb_freq, station0_, tile0_)
                               .Data());

      if (use_differential_beam_) {
        aocommon::MC2x2F::ATimesB(ant_buffer_ptr,
                                  inverse_central_gain_[job.buffer_offset],
                                  gain_matrix);
      } else {
        gain_matrix.AssignTo(ant_buffer_ptr);
      }
    }
  }
}
