#include "lofargrid.h"
#include "./../telescope/lofar.h"
#include "./../common/system.h"

#include <aocommon/imagecoordinates.h>
#include <cmath>
// #include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using namespace everybeam::gridded_response;

LOFARGrid::LOFARGrid(telescope::Telescope* telescope_ptr,
                     const coords::CoordinateSystem& coordinate_system)
    : GriddedResponse(telescope_ptr, coordinate_system) {
  const telescope::LOFAR& lofartelescope =
      static_cast<const telescope::LOFAR&>(*telescope_);
  size_t ncpus = common::System::ProcessorCount();

  // Set private members
  nthreads_ = std::min(ncpus, lofartelescope.nstations_);
  threads_.resize(nthreads_);
  delay_dir_ = lofartelescope.ms_properties_.delay_dir;
  tile_beam_dir_ = lofartelescope.ms_properties_.tile_beam_dir;
  subband_frequency_ = lofartelescope.ms_properties_.subband_freq;
  use_channel_frequency_ = lofartelescope.options_.useChannelFrequency;
};

bool LOFARGrid::CalculateStation(std::complex<float>* buffer, double time,
                                 double frequency, const size_t station_idx) {
  aocommon::Lane<Job> lane(nthreads_);
  lane_ = &lane;

  SetITRFVectors(time);
  double sb_freq = use_channel_frequency_ ? frequency : subband_frequency_;
  // Dummy calculation of gain matrix, needed for multi-threading
  telescope_->GetStation(station_idx)
      ->Response(time, frequency, station0_, sb_freq, station0_, tile0_);

  for (size_t i = 0; i != nthreads_; ++i) {
    threads_[i] =
        std::thread(&LOFARGrid::CalcThread, this, buffer, time, frequency);
  }

  for (size_t y = 0; y != height_; ++y) {
    lane.write(Job{.y = y, .antenna_idx = station_idx});
  }

  lane.write_end();
  for (size_t i = 0; i != nthreads_; ++i) threads_[i].join();
  return true;
}

bool LOFARGrid::CalculateAllStations(std::complex<float>* buffer, double time,
                                     double frequency) {
  aocommon::Lane<Job> lane(nthreads_);
  lane_ = &lane;

  SetITRFVectors(time);

  double sb_freq = use_channel_frequency_ ? frequency : subband_frequency_;

  // Dummy loop, needed for multi-threading
  for (size_t i = 0; i != telescope_->GetNrStations(); ++i) {
    telescope_->GetStation(i)->Response(time, frequency, station0_, sb_freq,
                                        station0_, tile0_);
  }

  // Prepare threads
  for (size_t i = 0; i != nthreads_; ++i) {
    threads_[i] =
        std::thread(&LOFARGrid::CalcThread, this, buffer, time, frequency);
  }

  for (size_t y = 0; y != height_; ++y) {
    for (size_t antenna_idx = 0; antenna_idx != telescope_->GetNrStations();
         ++antenna_idx) {
      lane.write(Job{.y = y, .antenna_idx = antenna_idx});
    }
  }

  lane.write_end();
  for (size_t i = 0; i != nthreads_; ++i) threads_[i].join();
  return true;
}

void LOFARGrid::SetITRFVectors(double time) {
  coords::ITRFConverter itrf_converter(time);
  coords::SetITRFVector(itrf_converter.toDirection(delay_dir_), station0_);
  coords::SetITRFVector(itrf_converter.toDirection(tile_beam_dir_), tile0_);

  const casacore::Unit rad_unit("rad");

  casacore::MDirection l_dir(
      casacore::MVDirection(casacore::Quantity(ra_ + M_PI / 2, rad_unit),
                            casacore::Quantity(0, rad_unit)),
      casacore::MDirection::J2000);
  coords::SetITRFVector(itrf_converter.toDirection(l_dir), l_vector_itrf_);

  casacore::MDirection m_dir(
      casacore::MVDirection(casacore::Quantity(ra_, rad_unit),
                            casacore::Quantity(dec_ + M_PI / 2, rad_unit)),
      casacore::MDirection::J2000);
  coords::SetITRFVector(itrf_converter.toDirection(m_dir), m_vector_itrf_);

  casacore::MDirection n_dir(
      casacore::MVDirection(casacore::Quantity(ra_, rad_unit),
                            casacore::Quantity(dec_, rad_unit)),
      casacore::MDirection::J2000);
  coords::SetITRFVector(itrf_converter.toDirection(n_dir), n_vector_itrf_);
}

void LOFARGrid::CalcThread(std::complex<float>* buffer, double time,
                           double frequency) {
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
          base_buffer + job.antenna_idx * values_per_ant;
      matrix22c_t gain_matrix = telescope_->GetStation(job.antenna_idx)
                                    ->Response(time, frequency, itrf_direction,
                                               sb_freq, station0_, tile0_);

      // (Optional) differential beam logic is handled at WSClean level
      ant_buffer_ptr[0] = gain_matrix[0][0];
      ant_buffer_ptr[1] = gain_matrix[0][1];
      ant_buffer_ptr[2] = gain_matrix[1][0];
      ant_buffer_ptr[3] = gain_matrix[1][1];
    }
  }
}