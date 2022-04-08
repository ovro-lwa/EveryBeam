// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "config.h"
#include "../aterms/klfittingaterm.h"
#include "../coords/coordutils.h"

#include <aocommon/fits/fitsreader.h>

using everybeam::aterms::KlFittingATerm;
using everybeam::coords::CoordinateSystem;

BOOST_AUTO_TEST_SUITE(testklfitting)

// This test is allowed to fail. It will be fixed in AST-852
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(verify, 1)
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

  for (size_t freq_idx = 0; freq_idx < n_frequencies; ++freq_idx) {
    for (size_t time_idx = 0; time_idx < n_times; ++time_idx) {
      size_t img_idx = time_idx * n_frequencies + freq_idx;
      for (size_t ant_idx = 0; ant_idx < n_antennas; ++ant_idx) {
        const size_t antenna_file_index = img_idx * n_antennas + ant_idx;
        for (size_t pol_idx = 0; pol_idx < n_corr; ++pol_idx) {
          size_t file_index = antenna_file_index * n_corr + pol_idx;
          std::complex<float>* antenna_buffer =
              aterm_expected.data() + file_index * single_fits_buffersize;
          reader.ReadIndex(antenna_buffer, file_index);
        }
      }
    }
  }

  // Create aterm with Karhunen-Loève base functions
  // In AST-852, verify that the minus sign in the DL/DL phase center is needed
  const CoordinateSystem coord_system = {
      screen_width, screen_height, ra, dec, dl, dm, -0.5 * dl, -0.5 * dm};
  KlFittingATerm kl_fitting_aterm(coord_system);
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
      // Since the Calculate function is not implemented yet, it will throw a
      // runtime error. It will be implemented in AST-852.
      BOOST_REQUIRE_THROW(
          kl_fitting_aterm.Calculate(
              kl_buffer, start_time + time_idx * time_step,
              start_frequency + freq_idx * frequency_step, 0, nullptr),
          std::runtime_error);
    }
  }
  // Verify that the aterms generated are the same
  // Due to huge size of the aterm buffer (17201280 elements), we only check the
  // match on a subset of it.
  const size_t positions_to_check = 100000;
  BOOST_CHECK_EQUAL_COLLECTIONS(
      kl_aterm.begin(), kl_aterm.begin() + positions_to_check,
      aterm_expected.begin(), aterm_expected.begin() + positions_to_check);
}

BOOST_AUTO_TEST_SUITE_END()