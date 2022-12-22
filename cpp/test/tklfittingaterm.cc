// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include <xtensor/xarray.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xnpy.hpp>

#include "config.h"
#include "../aterms/klfittingaterm.h"

#include <aocommon/fits/fitsreader.h>

using aocommon::CoordinateSystem;
using everybeam::aterms::KlFittingATerm;

BOOST_AUTO_TEST_SUITE(klfittingaterm)

BOOST_AUTO_TEST_CASE(verify, *boost::unit_test::tolerance(1e-8)) {
  // Create aterm from fits file

  aocommon::FitsReader reader(KL_SCREEN_PATH, true, true);

  const size_t screen_width = reader.ImageWidth();
  const size_t screen_height = reader.ImageHeight();
  const size_t n_frequencies = reader.NFrequencies();
  const size_t n_antennas = reader.NAntennas();
  const size_t n_times = reader.NTimesteps();
  const size_t n_corr = reader.NMatrixElements();
  const double ra = reader.PhaseCentreRA();
  const double dec = reader.PhaseCentreDec();
  const double dl = reader.PhaseCentreDL();
  const double dm = reader.PhaseCentreDM();

  const double start_time = reader.TimeDimensionStart();
  const double time_step = reader.TimeDimensionIncr();
  const double start_frequency = reader.FrequencyDimensionStart();
  const double frequency_step = reader.FrequencyDimensionIncr();

  const size_t aterm_buffersize = screen_width * screen_height * n_antennas *
                                  n_corr * n_times * n_frequencies;
  std::vector<std::complex<float>> aterm_expected(aterm_buffersize);

  // Read aterms from FITS file

  // Fits reader reads one screen at the time
  const unsigned int single_fits_buffersize = screen_width * screen_height;

  std::vector<float> antenna_buffer_real(single_fits_buffersize);
  std::vector<float> antenna_buffer_imag(single_fits_buffersize);

  for (size_t freq_idx = 0; freq_idx < n_frequencies; ++freq_idx) {
    for (size_t time_idx = 0; time_idx < n_times; ++time_idx) {
      size_t img_idx = time_idx * n_frequencies + freq_idx;
      for (size_t ant_idx = 0; ant_idx < n_antennas; ++ant_idx) {
        const size_t antenna_file_index = img_idx * n_antennas + ant_idx;
        for (size_t pol_idx = 0; pol_idx < 2; ++pol_idx) {
          size_t file_index = antenna_file_index * n_corr + pol_idx * 2;
          reader.ReadIndex(antenna_buffer_real.data(), file_index);
          file_index = antenna_file_index * n_corr + pol_idx * 2 + 1;
          reader.ReadIndex(antenna_buffer_imag.data(), file_index);
          for (size_t pixel_index = 0; pixel_index < single_fits_buffersize;
               ++pixel_index) {
            // pol_idx is either 0 (for the XX element)
            // or 1 (for the YY element). These elements are placed on the
            // diagonal of a Jones matrix (2x2 complex matrix = 4 complex
            // elements). When using a single index for the Jones matrix the
            // index of the first diagonal element is 0, and of the
            // second diagonal element is 3.
            constexpr int kNrElementsPerJonesMatrix = 4;
            constexpr int kDistanceBetweenDiagonalElements = 3;
            aterm_expected[antenna_file_index * n_corr *
                               single_fits_buffersize +
                           pol_idx * kDistanceBetweenDiagonalElements +
                           kNrElementsPerJonesMatrix * pixel_index] =
                std::complex<float>(antenna_buffer_real[pixel_index],
                                    antenna_buffer_imag[pixel_index]);
          }
        }
      }
    }
  }

  // Create aterm with Karhunen-Loève base functions
  const CoordinateSystem coord_system = {
      screen_width, screen_height, ra, dec, dl, dm, 0.0, 0.0};
  // Empty selection of stations will select all stations in h5parm
  std::vector<std::string> station_names = {};
  const int kOrder = 6;
  KlFittingATerm kl_fitting_aterm(station_names, coord_system, kOrder);
  kl_fitting_aterm.Open({H5_SOLUTIONS_PATH});

  std::vector<std::complex<float>> kl_aterm(aterm_buffersize);

  // KL aterms generator will produce aterms for a given time and frequency.
  // It handles multiple antennas and polarizations internally.
  const unsigned int single_kl_buffersize =
      screen_width * screen_height * n_corr * n_antennas;

  for (size_t freq_idx = 0; freq_idx < n_frequencies; ++freq_idx) {
    for (size_t time_idx = 0; time_idx < n_times; ++time_idx) {
      size_t img_idx = time_idx * n_frequencies + freq_idx;
      std::complex<float>* kl_buffer =
          kl_aterm.data() + img_idx * single_kl_buffersize;
      kl_fitting_aterm.Calculate(kl_buffer, start_time + time_idx * time_step,
                                 start_frequency + freq_idx * frequency_step, 0,
                                 nullptr);
    }
  }
  // Verify that the aterms generated are the same
  const size_t positions_to_check = kl_aterm.size();
  for (size_t i = 0; i < positions_to_check; ++i) {
    // Compare the aterm computed by the KlFittingATerm (C++ class)
    // to the one computed by the python script from the ska-sdp-screen-fitting
    // repository Note that the python script uses a slightly different method
    // for interpolation. That method is prone to phase wrapping errors, so we
    // do not want to use it for KlFittingAterm. Due to the difference in
    // methods the threshold needs to be higher than usual
    BOOST_CHECK(abs(kl_aterm[i] - aterm_expected[i]) < 0.2);
  }
}

BOOST_AUTO_TEST_SUITE_END()
