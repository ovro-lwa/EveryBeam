// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "aterms/atermconfig.h"

#include <algorithm>
#include <map>
#include <string>

#include <aocommon/fits/fitsreader.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <schaapcommon/h5parm/h5parm.h>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/tokenizer.hpp>

#include "config.h"
#include "load.h"
#include "aterms/parsetprovider.h"

using casacore::ArrayColumn;
using casacore::ArrayMeasColumn;
using casacore::IPosition;
using casacore::MDirection;
using casacore::Table;

using everybeam::aterms::ATermConfig;

const size_t kAtermMatrixSize = 4;

namespace {

class MockParsetProvider : public everybeam::aterms::ParsetProvider {
 public:
  MockParsetProvider(bool provide_beam, bool provide_diagonal,
                     bool provide_fourier_fit, bool provide_kl_fit)
      : provide_beam_(provide_beam),
        provide_diagonal_(provide_diagonal),
        provide_fourier_fit_(provide_fourier_fit),
        provide_kl_fit_(provide_kl_fit) {}

  std::string GetString(const std::string& key) const final override {
    BOOST_FAIL("GetString: Unexpected key: " + key);
    return "";
  };

  std::string GetStringOr([[maybe_unused]] const std::string& key,
                          const std::string& orValue) const final override {
    if (provide_fourier_fit_ || provide_kl_fit_) {
      if ((key == "fourierfit.solutions") || (key == "klfit.solutions")) {
        return H5_SOLUTIONS_PATH;
      }
    }
    return orValue;
  };

  std::vector<std::string> GetStringList(
      const std::string& key) const final override {
    std::vector<std::string> list;

    if (key == "aterms") {
      if (provide_beam_) {
        list.push_back("beam");
      }
      if (provide_diagonal_) {
        list.push_back("diagonal");
      }
      if (provide_fourier_fit_) {
        list.push_back("fourierfit");
      }
      if (provide_kl_fit_) {
        list.push_back("klfit");
      }
      return list;
    }

    if (provide_diagonal_) {
      if (key == "diagonal.images") {
        list.push_back(KL_SCREEN_PATH);
        return list;
      }
    }

    BOOST_FAIL("GetStringList: Unexpected key: " + key);
    return list;
  };

  bool GetBool(const std::string& key) const final override {
    BOOST_FAIL("GetBool: Unexpected key: " + key);
    return false;
  }

  bool GetBoolOr(const std::string& key, bool orValue) const final override {
    if (provide_beam_) {
      if (key == "beam.differential") {
        return true;
      } else if (key == "beam.usechannelfreq") {
        return true;
      } else if (key == "beam.frequency_interpolation") {
        return true;
      }
    }

    return orValue;
  }

  double GetDoubleOr(const std::string& key,
                     double orValue) const final override {
    if (provide_beam_) {
      if (key == "beam.update_interval") {
        return 1.0;
      }
    }

    return orValue;
  };

 private:
  const bool provide_beam_;
  const bool provide_diagonal_;
  const bool provide_fourier_fit_;
  const bool provide_kl_fit_;
};

}  // namespace

void GenerateAterms(const casacore::MeasurementSet& ms,
                    const MockParsetProvider& parset_provider,
                    std::vector<std::complex<float>>& aterm_buffer,
                    double start_time, double time_step, double start_frequency,
                    double frequency_step, size_t n_antennas, size_t n_channels,
                    size_t n_times,
                    everybeam::coords::CoordinateSystem system) {
  everybeam::ATermSettings aterm_settings;
  bool is_updated = false;
  ATermConfig aterm(n_antennas, system, aterm_settings);
  aterm.SetSaveATerms(false, "");
  aterm.Read(ms, parset_provider, "");

  const size_t single_buffer_size =
      system.width * system.height * n_antennas * kAtermMatrixSize;
  size_t prev_image_index = 0;
  // Frequency / timesteps are from the MS
  for (size_t freq_index = 0; freq_index < n_channels; ++freq_index) {
    for (size_t time_index = 0; time_index < n_times; ++time_index) {
      size_t image_index = time_index * n_channels + freq_index;
      std::complex<float>* buffer =
          aterm_buffer.data() + image_index * single_buffer_size;
      is_updated = aterm.Calculate(
          buffer, start_time + time_index * time_step,
          start_frequency + freq_index * frequency_step, 0, nullptr);
      if (!is_updated) {
        BOOST_CHECK(time_index > 0);
        prev_image_index = (time_index - 1) * n_channels + freq_index;
        std::copy_n(aterm_buffer.data() + prev_image_index * single_buffer_size,
                    single_buffer_size,
                    aterm_buffer.data() + image_index * single_buffer_size);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE(aterm_config)

BOOST_AUTO_TEST_CASE(combine_aterms) {
  casacore::MeasurementSet ms(SCREEN_FITTING_MS);

  // Get data from ms
  Table antenna_table(ms.keywordSet().asTable("ANTENNA"));
  const size_t n_antennas_ms = antenna_table.nrow();
  Table field_table(ms.keywordSet().asTable("FIELD"));
  ArrayMeasColumn<MDirection> field_column_1(field_table, "PHASE_DIR");
  const MDirection phase_center_ms = *(field_column_1(0).data());

  // Get data from fits file
  aocommon::FitsReader fitsreader(KL_SCREEN_PATH, true, true);
  const double frequency_step = fitsreader.FrequencyDimensionIncr();
  const float start_frequency = fitsreader.FrequencyDimensionStart();
  const size_t n_channels = fitsreader.NFrequencies();
  const size_t n_times = fitsreader.NTimesteps();
  const float start_time = fitsreader.TimeDimensionStart();
  const double time_step = fitsreader.TimeDimensionIncr();
  const size_t n_antennas = fitsreader.NAntennas();
  const double ra = fitsreader.PhaseCentreRA();
  const double dec = fitsreader.PhaseCentreDec();
  const double dl = fitsreader.PhaseCentreDL();
  const double dm = fitsreader.PhaseCentreDM();
  const size_t screen_width = fitsreader.ImageWidth();
  const size_t screen_height = fitsreader.ImageHeight();

  // Verify that the MS and the fits file contain the same information
  BOOST_TEST(n_antennas_ms == n_antennas);
  BOOST_TEST(phase_center_ms.getAngle("rad").getValue()[0] == ra,
             boost::test_tools::tolerance(0.001));
  BOOST_TEST(phase_center_ms.getAngle("rad").getValue()[1] == dec,
             boost::test_tools::tolerance(0.001));

  // Initialize coordinate system
  const everybeam::coords::CoordinateSystem coord_system = {
      screen_width, screen_height, ra, dec, dl, dm, -0.5 * dl, -0.5 * dm};

  // Create parset providers
  const MockParsetProvider provider_combined(true, true, false, false);
  const MockParsetProvider provider_beam(true, false, false, false);
  const MockParsetProvider provider_diagonal(false, true, false, false);

  // Initialize aterm buffers
  const size_t aterm_buffer_size = coord_system.width * coord_system.height *
                                   n_antennas * kAtermMatrixSize * n_channels *
                                   n_times;
  std::vector<std::complex<float>> aterm_buffer_combined(aterm_buffer_size);
  std::vector<std::complex<float>> aterm_buffer_beam(aterm_buffer_size);
  std::vector<std::complex<float>> aterm_buffer_diagonal(aterm_buffer_size);

  // Generate aterms for the three different configurations
  GenerateAterms(ms, provider_combined, aterm_buffer_combined, start_time,
                 time_step, start_frequency, frequency_step, n_antennas,
                 n_channels, n_times, coord_system);
  GenerateAterms(ms, provider_beam, aterm_buffer_beam, start_time, time_step,
                 start_frequency, frequency_step, n_antennas, n_channels,
                 n_times, coord_system);
  GenerateAterms(ms, provider_diagonal, aterm_buffer_diagonal, start_time,
                 time_step, start_frequency, frequency_step, n_antennas,
                 n_channels, n_times, coord_system);

  // Combine the aterms generated in the beam-only and diagonal-only modes
  std::vector<std::complex<float>> aterms_product(aterm_buffer_size);
  for (size_t j = 0; j != aterm_buffer_size; j += kAtermMatrixSize) {
    aocommon::MC2x2F scratch;
    aocommon::Matrix2x2::ATimesB(scratch.Data(),
                                 aterm_buffer_diagonal.data() + j,
                                 aterm_buffer_beam.data() + j);
    aocommon::Matrix2x2::Assign(aterms_product.data() + j, scratch.Data());
  }

  BOOST_CHECK_EQUAL_COLLECTIONS(aterms_product.begin(), aterms_product.end(),
                                aterm_buffer_combined.begin(),
                                aterm_buffer_combined.end());
}

BOOST_AUTO_TEST_CASE(fourierfit) {
  casacore::MeasurementSet ms(SCREEN_FITTING_MS);

  // Get data from ms
  Table antenna_table(ms.keywordSet().asTable("ANTENNA"));
  const size_t n_antennas_ms = antenna_table.nrow();
  Table field_table(ms.keywordSet().asTable("FIELD"));
  ArrayMeasColumn<MDirection> field_column_1(field_table, "PHASE_DIR");
  const MDirection phase_center_ms = *(field_column_1(0).data());
  const double ra = phase_center_ms.getAngle("rad").getValue()[0];
  const double dec = phase_center_ms.getAngle("rad").getValue()[1];

  // Get data from h5parm file
  schaapcommon::h5parm::H5Parm h5parmfile(H5_SOLUTIONS_PATH);
  schaapcommon::h5parm::SolTab phase_soltab = h5parmfile.GetSolTab("phase000");
  std::vector<double> frequencies = phase_soltab.GetRealAxis("freq");
  std::vector<double> times = phase_soltab.GetRealAxis("time");
  double frequency_step = phase_soltab.GetFreqInterval();
  double time_step = phase_soltab.GetTimeInterval();
  const size_t n_channels = frequencies.size();
  const size_t n_times = times.size();
  unsigned int n_antennas = phase_soltab.GetAxis("ant").size;

  // Verify that the MS and the h5 file contain the same information
  BOOST_TEST(n_antennas_ms == n_antennas);

  // Properties of grid
  const size_t width = 32;
  const size_t height = 32;
  const double dl = 3.5 * M_PI / 180. / width;
  const double dm = 3.5 * M_PI / 180. / height;

  // Initialize coordinate system
  const everybeam::coords::CoordinateSystem coord_system = {
      width, height, ra, dec, dl, dm, 0., 0.};

  // Create parset providers
  const MockParsetProvider provider_fourierfit(false, false, true, false);

  // Initialize aterm buffers
  const size_t aterm_buffer_size = coord_system.width * coord_system.height *
                                   n_antennas * kAtermMatrixSize * n_channels *
                                   n_times;
  std::vector<std::complex<float>> aterm_buffer_fourierfit(aterm_buffer_size);

  GenerateAterms(ms, provider_fourierfit, aterm_buffer_fourierfit, times[0],
                 time_step, frequencies[0], frequency_step, n_antennas,
                 n_channels, n_times, coord_system);
}

BOOST_AUTO_TEST_SUITE_END()