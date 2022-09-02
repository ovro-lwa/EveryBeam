// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "phasedarraygrid.h"
#include "../telescope/phasedarray.h"

#include <aocommon/lane.h>
#include <aocommon/imagecoordinates.h>
#include <aocommon/threadpool.h>
#include <cmath>
#include <iostream>
namespace everybeam {

using telescope::PhasedArray;

namespace griddedresponse {

PhasedArrayGrid::PhasedArrayGrid(
    const telescope::Telescope* telescope_ptr,
    const coords::CoordinateSystem& coordinate_system)
    : GriddedResponse(telescope_ptr, coordinate_system),
      PhasedArrayResponse(static_cast<const PhasedArray*>(telescope_ptr)) {
  // Compute and set number of threads
  const size_t ncpus = aocommon::system::ProcessorCount();
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
  const telescope::PhasedArray& phased_array =
      static_cast<const telescope::PhasedArray&>(*telescope_);
  aocommon::Lane<Job> lane(threads_.size());
  lane_ = &lane;

  SetITRFVectors(time);

  bool apply_normalisation = false;
  inverse_central_gain_.resize(phased_array.GetNrStations());
  for (size_t i = 0; i != phased_array.GetNrStations(); ++i) {
    apply_normalisation = CalculateBeamNormalisation(
        beam_mode, time, frequency, i, inverse_central_gain_[i]);
  }

  // Prepare threads
  for (auto& thread : threads_) {
    thread = std::thread(&PhasedArrayGrid::CalcThread, this, beam_mode,
                         apply_normalisation, buffer, time, frequency);
  }

  for (size_t y = 0; y != height_; ++y) {
    for (size_t antenna_idx = 0; antenna_idx != phased_array.GetNrStations();
         ++antenna_idx) {
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

void PhasedArrayGrid::CalcThread(BeamMode beam_mode, bool apply_normalisation,
                                 std::complex<float>* buffer, double time,
                                 double frequency) {
  const telescope::PhasedArray& phased_array =
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
      const double sqrt_term = 1.0 - l * l - m * m;
      if (sqrt_term >= 0.0) {
        n = std::sqrt(sqrt_term);
      } else {
        n = -std::sqrt(-sqrt_term);
      }

      const vector3r_t itrf_direction =
          l * l_vector_itrf_ + m * m_vector_itrf_ + n * n_vector_itrf_;

      std::complex<float>* base_buffer = buffer + (x + job.y * width_) * 4;

      std::complex<float>* ant_buffer_ptr =
          base_buffer + job.buffer_offset * values_per_ant;

      const aocommon::MC2x2F gain_matrix = aocommon::MC2x2F(
          phased_array.GetStation(job.antenna_idx)
              .Response(beam_mode, time, frequency, itrf_direction, sb_freq,
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
}  // namespace griddedresponse
}  // namespace everybeam
