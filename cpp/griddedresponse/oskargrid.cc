#include "oskargrid.h"
#include "./../telescope/oskar.h"

#include <aocommon/imagecoordinates.h>
#include <aocommon/threadpool.h>
#include <cmath>
#include <iostream>

using everybeam::griddedresponse::OSKARGrid;

OSKARGrid::OSKARGrid(telescope::Telescope* telescope_ptr,
                     const coords::CoordinateSystem& coordinate_system)
    : GriddedResponse(telescope_ptr, coordinate_system) {
  const telescope::OSKAR& oskartelescope =
      dynamic_cast<const telescope::OSKAR&>(*telescope_);
  size_t ncpus = aocommon::ThreadPool::NCPUs();

  // Set private members
  nthreads_ = std::min(ncpus, oskartelescope.nstations_);
  threads_.resize(nthreads_);
  delay_dir_ = oskartelescope.ms_properties_.delay_dir;
  tile_beam_dir_ = oskartelescope.ms_properties_.delay_dir;
};

void OSKARGrid::CalculateStation(std::complex<float>* buffer, double time,
                                 double frequency, size_t station_idx, size_t) {
  const telescope::OSKAR& oskartelescope =
      static_cast<const telescope::OSKAR&>(*telescope_);
  aocommon::Lane<Job> lane(nthreads_);
  lane_ = &lane;

  SetITRFVectors(time);

  for (size_t i = 0; i != nthreads_; ++i) {
    threads_[i] =
        std::thread(&OSKARGrid::CalcThread, this, buffer, time, frequency);
  }

  for (size_t y = 0; y != height_; ++y) {
    lane.write(Job{.y = y, .antenna_idx = station_idx, .buffer_offset = 0});
  }

  lane.write_end();
  for (size_t i = 0; i != nthreads_; ++i) threads_[i].join();
}

void OSKARGrid::CalculateAllStations(std::complex<float>* buffer, double time,
                                     double frequency, size_t) {
  const telescope::OSKAR& oskartelescope =
      static_cast<const telescope::OSKAR&>(*telescope_);
  aocommon::Lane<Job> lane(nthreads_);
  lane_ = &lane;

  SetITRFVectors(time);

  // Prepare threads
  for (size_t i = 0; i != nthreads_; ++i) {
    threads_[i] =
        std::thread(&OSKARGrid::CalcThread, this, buffer, time, frequency);
  }

  for (size_t y = 0; y != height_; ++y) {
    for (size_t antenna_idx = 0; antenna_idx != oskartelescope.GetNrStations();
         ++antenna_idx) {
      lane.write(Job{
          .y = y, .antenna_idx = antenna_idx, .buffer_offset = antenna_idx});
    }
  }

  lane.write_end();
  for (size_t i = 0; i != nthreads_; ++i) threads_[i].join();
}

void OSKARGrid::SetITRFVectors(double time) {
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

void OSKARGrid::CalcThread(std::complex<float>* buffer, double time,
                           double frequency) {
  const telescope::OSKAR& oskartelescope =
      static_cast<const telescope::OSKAR&>(*telescope_);
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

      std::complex<float>* base_buffer = buffer + (x + job.y * height_) * 4;

      std::complex<float>* ant_buffer_ptr =
          base_buffer + job.buffer_offset * values_per_ant;

      matrix22c_t gain_matrix = oskartelescope.GetStation(job.antenna_idx)
                                    ->Response(time, frequency, itrf_direction,
                                               sb_freq, station0_, tile0_);
      ant_buffer_ptr[0] = gain_matrix[0][0];
      ant_buffer_ptr[1] = gain_matrix[0][1];
      ant_buffer_ptr[2] = gain_matrix[1][0];
      ant_buffer_ptr[3] = gain_matrix[1][1];
    }
  }
}
