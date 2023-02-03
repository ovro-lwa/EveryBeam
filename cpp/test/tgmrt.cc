// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "config.h"

#include "../load.h"
#include "../griddedresponse/dishgrid.h"
#include "../telescope/telescope.h"
#include "../telescope/dish.h"

namespace {

const double kTime = 5.08859e+09;
// There's a "reference" fits file at 383.325 MHz that I've used to
// determine some properties such as the (approx) FWHM that are later
// used in these tests
const double kFrequency = 383.325e6;
const size_t kWidth = 64;
const size_t kHeight = 64;

aocommon::CoordinateSystem GetCoordinateSystem() {
  aocommon::CoordinateSystem coord_system;
  coord_system.width = kWidth;
  coord_system.height = kHeight;
  coord_system.ra = -156.837226063 * (M_PI / 180.0);
  coord_system.dec = 50.4691390331 * (M_PI / 180.0);
  coord_system.dl = 1.0 * (M_PI / 180.0 / 60.0);
  coord_system.dm = 1.0 * (M_PI / 180.0 / 60.0);
  coord_system.l_shift = 0.0;
  coord_system.m_shift = 0.0;
  return coord_system;
}

struct GmrtFixture {
  GmrtFixture() : ms(GMRT_MOCK_MS), telescope(everybeam::Load(ms, options)) {
    grid_response = telescope->GetGriddedResponse(GetCoordinateSystem());
  }
  ~GmrtFixture(){};
  casacore::MeasurementSet ms;
  everybeam::Options options;
  std::unique_ptr<everybeam::telescope::Telescope> telescope;
  std::unique_ptr<everybeam::griddedresponse::GriddedResponse> grid_response;

  size_t Index(size_t x, size_t y) const { return 4 * (kWidth * y + x); }

  void CheckStation(const std::complex<float>* single_antenna_buffer) const {
    const size_t size = kWidth * kHeight * 2 * 2;
    for (size_t i = 0; i != size; i += 4) {
      // Check if off diagonals are zero:
      BOOST_CHECK_LT(std::fabs(single_antenna_buffer[i + 1].real()), 1e-4);
      BOOST_CHECK_LT(std::fabs(single_antenna_buffer[i + 2].real()), 1e-4);
      // Because X and Y are calculated in the same way, they should always be
      // identical
      BOOST_CHECK_CLOSE_FRACTION(single_antenna_buffer[i].real(),
                                 single_antenna_buffer[i + 3].real(), 1e-4);
    }
    // Check if all phases (/imaginary values) are zero:
    for (size_t i = 0; i != size; ++i) {
      BOOST_CHECK_LT(std::fabs(single_antenna_buffer[i].imag()), 1e-4);
    }

    // Check if the centre pixel is approximately the unit matrix
    // (Due to interpolation the diagonal gains are not exactly 1; it was
    // 0.999789 when I checked).
    const size_t centre_index = Index(kWidth / 2, kHeight / 2);
    BOOST_CHECK_CLOSE_FRACTION(single_antenna_buffer[centre_index].real(), 1.0,
                               1e-2);

    // The FWHM is approx at the edge. These are (un-squared) voltages, hence
    // its value is about sqrt(0.5) (=~0.7).
    const float value_at_edge = 0.697981417;
    const size_t edge_index = Index(kWidth / 2, 0);
    BOOST_CHECK_CLOSE_FRACTION(single_antenna_buffer[edge_index].real(),
                               value_at_edge, 1e-2);
    BOOST_CHECK_CLOSE_FRACTION(single_antenna_buffer[edge_index + 3].real(),
                               value_at_edge, 1e-2);
  }
};

void CheckUndersampledIntegratedResponse(size_t undersampling_factor) {
  GmrtFixture gmrt;
  const std::vector<double> time_array{kTime};
  const size_t n_antennas = gmrt.telescope->GetNrStations();
  const size_t n_baselines = n_antennas * (n_antennas + 1) / 2;
  const std::vector<double> baseline_weights(n_baselines * time_array.size(),
                                             1.0);

  const std::vector<aocommon::HMC4x4> result =
      gmrt.grid_response->UndersampledIntegratedResponse(
          everybeam::CorrectionMode::kFull, time_array, kFrequency, 0,
          undersampling_factor, baseline_weights);
  const size_t undersampled_width = kWidth / undersampling_factor;
  const size_t undersampled_height = kHeight / undersampling_factor;
  BOOST_REQUIRE_EQUAL(result.size(), undersampled_width * undersampled_height);

  for (size_t element_index = 0; element_index != 4; ++element_index) {
    const float lowres_centre_pixel =
        result[undersampled_width * (undersampled_height / 2) +
               undersampled_width / 2]
            .Data(element_index);
    const size_t off_centre_y =
        undersampled_height / 2 - 8 / undersampling_factor;
    const float lowres_off_centre_pixel =
        result[undersampled_width * off_centre_y + undersampled_width / 2].Data(
            element_index);
    const float lowres_edge_pixel =
        result[undersampled_width / 2].Data(element_index);

    std::vector<float> buffer(kWidth * kHeight);
    gmrt.grid_response->UpsampleResponse(buffer.data(), element_index, kWidth,
                                         kHeight, result, undersampling_factor);

    const float upsampled_centre_pixel =
        buffer[kWidth * (kHeight / 2) + kWidth / 2];
    const float upsampled_edge_pixel = buffer[kWidth / 2];

    if (element_index == 0 || element_index == 3) {
      BOOST_CHECK_CLOSE_FRACTION(lowres_centre_pixel, 1.0, 1e-2);
      BOOST_CHECK_CLOSE_FRACTION(lowres_off_centre_pixel, 0.959, 1e-2);
      BOOST_CHECK_CLOSE_FRACTION(lowres_edge_pixel, 0.487, 1e-2);
      BOOST_CHECK_CLOSE_FRACTION(upsampled_centre_pixel, 1.0, 1e-2);
      BOOST_CHECK_CLOSE_FRACTION(upsampled_edge_pixel, 0.487, 1e-2);
    } else {
      BOOST_CHECK_LT(std::fabs(lowres_centre_pixel), 1e-2);
      BOOST_CHECK_LT(std::fabs(lowres_off_centre_pixel), 1e-2);
      BOOST_CHECK_LT(std::fabs(lowres_edge_pixel), 1e-2);
      BOOST_CHECK_LT(std::fabs(upsampled_centre_pixel), 1e-2);
      BOOST_CHECK_LT(std::fabs(upsampled_edge_pixel), 1e-2);
    }
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(tgmrt)

BOOST_FIXTURE_TEST_CASE(single_antenna_gridded_response, GmrtFixture) {
  BOOST_CHECK(nullptr != dynamic_cast<everybeam::griddedresponse::DishGrid*>(
                             grid_response.get()));

  // Define buffer and get gridded responses
  std::vector<std::complex<float>> single_antenna_buffer(
      grid_response->GetStationBufferSize(1));
  grid_response->Response(everybeam::CorrectionMode::kFull,
                          single_antenna_buffer.data(), kTime, kFrequency, 0,
                          0);
  BOOST_REQUIRE_EQUAL(single_antenna_buffer.size(),
                      std::size_t(kWidth * kHeight * 2 * 2));

  CheckStation(single_antenna_buffer.data());
}

BOOST_FIXTURE_TEST_CASE(all_antennas_gridded_response, GmrtFixture) {
  std::vector<std::complex<float>> buffer(
      grid_response->GetStationBufferSize(telescope->GetNrStations()));
  grid_response->ResponseAllStations(everybeam::CorrectionMode::kFull,
                                     buffer.data(), kTime, kFrequency, 0);
  const size_t single_antenna_size = kWidth * kHeight * 4;
  const size_t expected_size = telescope->GetNrStations() * single_antenna_size;
  BOOST_REQUIRE_EQUAL(buffer.size(), expected_size);

  for (size_t antenna_index = 0; antenna_index != telescope->GetNrStations();
       ++antenna_index) {
    const std::complex<float>* single_antenna_buffer =
        &buffer[single_antenna_size * antenna_index];
    CheckStation(single_antenna_buffer);
  }
}

BOOST_AUTO_TEST_CASE(undersampled_integrated_response_factor_1) {
  CheckUndersampledIntegratedResponse(1);
}

BOOST_AUTO_TEST_CASE(undersampled_integrated_response_factor_2) {
  CheckUndersampledIntegratedResponse(2);
}

BOOST_AUTO_TEST_SUITE_END()
