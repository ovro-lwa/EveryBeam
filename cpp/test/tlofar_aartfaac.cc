// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "config.h"
#include "../beammode.h"
#include "../telescope/lofar.h"
#include "../griddedresponse/griddedresponse.h"
#include "../griddedresponse/aartfaacgrid.h"
#include "../pointresponse/aartfaacpoint.h"
#include "../elementresponse.h"
#include "../load.h"

using everybeam::BeamMode;
using everybeam::ElementResponseModel;
using everybeam::Load;
using everybeam::Options;
using everybeam::coords::CoordinateSystem;
using everybeam::griddedresponse::AartfaacGrid;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::pointresponse::AartfaacPoint;
using everybeam::telescope::LOFAR;
using everybeam::telescope::Telescope;

BOOST_AUTO_TEST_SUITE(tlofar_aartfaac)

namespace {
// Properties extracted from MS
const double kTime = 5.1039942e+09;
const double kFrequency = 32617187.5;
const double kRa = -0.50402002;
const double kDec = 0.91858207;

// Properties of grid
const size_t kWidth = 4;
const size_t kHeight = 4;
const double kDl = 0.5 * M_PI / 180.;
const double kDm = 0.5 * M_PI / 180.;
const double kShiftL = 0.0;
const double kShiftM = 0.0;
const CoordinateSystem kCoordSystem = {kWidth, kHeight, kRa,     kDec,
                                       kDl,    kDm,     kShiftL, kShiftM};

std::unique_ptr<Telescope> LoadAartfaacTelescope(const Options& options) {
  std::unique_ptr<Telescope> telescope = Load(AARTFAAC_LBA_MOCK_MS, options);
  const size_t n_stations = telescope->GetNrStations();
  BOOST_CHECK_EQUAL(n_stations, 576);
  return telescope;
}
}  // namespace

BOOST_AUTO_TEST_CASE(load) {
  Options options;
  {
    options.element_response_model = ElementResponseModel::kLOBES;
    BOOST_CHECK_THROW(LoadAartfaacTelescope(options), std::runtime_error);
  }
  {
    options.element_response_model = ElementResponseModel::kHamakerLba;
    auto telescope = LoadAartfaacTelescope(options);
    BOOST_CHECK(dynamic_cast<LOFAR*>(telescope.get()));
  }
  {
    options.element_response_model = ElementResponseModel::kHamakerLba;
    auto telescope = LoadAartfaacTelescope(options);
    BOOST_CHECK(dynamic_cast<LOFAR*>(telescope.get()));
  }
}

BOOST_AUTO_TEST_CASE(gridded_response, *boost::unit_test::tolerance(1e-8)) {
  const Options options;
  auto telescope = LoadAartfaacTelescope(options);

  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(kCoordSystem);
  BOOST_CHECK(dynamic_cast<AartfaacGrid*>(grid_response.get()));

  const size_t buffer_size = grid_response->GetStationBufferSize(1);
  BOOST_CHECK_EQUAL(buffer_size, 4u * kWidth * kHeight);

  std::vector<std::complex<float>> full_response(buffer_size);
  std::vector<std::complex<float>> array_factor(buffer_size);
  std::vector<std::complex<float>> element_response(buffer_size);

  // Reference solution element response
  std::vector<std::complex<float>> element_response_ref(buffer_size);
  grid_response->Response(everybeam::BeamMode::kElement,
                          element_response_ref.data(), kTime, kFrequency, 0, 0);

  // Reference solution array factor (identity)
  const std::vector<std::complex<float>> array_factor_ref = {
      {1, 0}, {0, 0}, {0, 0}, {1, 0}};

  // Compute station responses all at once
  const size_t buffer_size_all =
      grid_response->GetStationBufferSize(telescope->GetNrStations());
  std::vector<std::complex<float>> full_response_all(buffer_size_all);
  grid_response->ResponseAllStations(BeamMode::kFull, full_response_all.data(),
                                     kTime, kFrequency, 0);

  for (size_t i = 0; i < telescope->GetNrStations(); ++i) {
    grid_response->Response(BeamMode::kFull, full_response.data(), kTime,
                            kFrequency, i, 0);
    grid_response->Response(BeamMode::kArrayFactor, array_factor.data(), kTime,
                            kFrequency, i, 0);
    grid_response->Response(BeamMode::kElement, element_response.data(), kTime,
                            kFrequency, i, 0);

    // Element response should be identical to reference element response
    BOOST_CHECK_EQUAL_COLLECTIONS(
        element_response.data(), element_response.data() + buffer_size,
        element_response_ref.data(), element_response_ref.data() + buffer_size);

    // Full response should be identical to element response
    BOOST_CHECK_EQUAL_COLLECTIONS(
        full_response.data(), full_response.data() + buffer_size,
        element_response.data(), element_response.data() + buffer_size);

    // Full response for single station should match corresponding part
    // in response for all stations
    BOOST_CHECK_EQUAL_COLLECTIONS(
        full_response_all.data() + i * buffer_size,
        full_response_all.data() + (i + 1) * buffer_size, full_response.data(),
        full_response.data() + buffer_size);

    for (size_t pixel = 0; pixel < kWidth * kHeight; ++pixel)
      // Array factor should be identity
      BOOST_CHECK_EQUAL_COLLECTIONS(array_factor.data() + pixel * 4u,
                                    array_factor.data() + (pixel + 1u) * 4u,
                                    array_factor_ref.data(),
                                    array_factor_ref.data() + 4);
  }
}

BOOST_AUTO_TEST_CASE(point_response, *boost::unit_test::tolerance(1e-8)) {
  const Options options;
  auto telescope = LoadAartfaacTelescope(options);

  const size_t station_index = 7;
  // Compute station response for station 7
  auto grid_response = telescope->GetGriddedResponse(kCoordSystem);
  const size_t buffer_size = grid_response->GetStationBufferSize(1);
  std::vector<std::complex<float>> full_response(buffer_size);
  grid_response->Response(BeamMode::kFull, full_response.data(), kTime,
                          kFrequency, station_index, 0);

  // Compute point response
  auto point_response = telescope->GetPointResponse(kTime);
  BOOST_CHECK(dynamic_cast<AartfaacPoint*>(point_response.get()));

  std::vector<std::complex<float>> point_buffer_all_stations(
      point_response->GetAllStationsBufferSize());
  point_response->ResponseAllStations(
      BeamMode::kFull, point_buffer_all_stations.data(), kCoordSystem.ra,
      kCoordSystem.dec, kFrequency, 0);

  // Offset for center pixel
  std::size_t offset_22 = (2 + 2 * kWidth) * 4;
  // Offset for station 7 in point response
  std::size_t offset_point = 4u * station_index;

  BOOST_CHECK_EQUAL_COLLECTIONS(
      point_buffer_all_stations.data() + offset_point,
      point_buffer_all_stations.data() + offset_point + 4,
      full_response.data() + offset_22, full_response.data() + offset_22 + 4);
}

BOOST_AUTO_TEST_SUITE_END()