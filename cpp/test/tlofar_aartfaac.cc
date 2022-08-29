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
using everybeam::pointresponse::PointResponse;
using everybeam::telescope::LOFAR;
using everybeam::telescope::Telescope;

BOOST_AUTO_TEST_SUITE(tlofar_aartfaac)

namespace {
// Properties extracted from MS
const double kTime = 5.1039942e+09;
const double kFrequency = 32617187.5;
const double kRa = -0.50402002;
const double kDec = 0.91858207;

const everybeam::vector3r_t kDirection{0.663096, -0.0590573, 0.746199};

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

static std::unique_ptr<Telescope> LoadAartfaac6Telescope(
    const Options& options) {
  std::unique_ptr<Telescope> telescope = Load(AARTFAAC_6_LBA_MOCK_MS, options);
  const size_t n_stations = telescope->GetNrStations();
  BOOST_REQUIRE_EQUAL(n_stations, 288);
  return telescope;
}

BOOST_AUTO_TEST_CASE(load) {
  Options options;
  // TODO AST-807 harden this test
  // - When an Aartfaac-12 MS is used it throws.
  {
    options.element_response_model = ElementResponseModel::kLOBES;
    std::unique_ptr<Telescope> telescope = LoadAartfaac6Telescope(options);
    BOOST_CHECK(dynamic_cast<LOFAR*>(telescope.get()));
  }
  {
    options.element_response_model = ElementResponseModel::kHamakerLba;
    std::unique_ptr<Telescope> telescope = LoadAartfaacTelescope(options);
    BOOST_CHECK(dynamic_cast<LOFAR*>(telescope.get()));
  }
}

BOOST_AUTO_TEST_CASE(gridded_response_hamaker,
                     *boost::unit_test::tolerance(1e-8)) {
  Options options;
  options.element_response_model = ElementResponseModel::kHamakerLba;
  std::unique_ptr<Telescope> telescope = LoadAartfaacTelescope(options);

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
  std::unique_ptr<Telescope> telescope = LoadAartfaacTelescope(options);

  const size_t station_index = 7;
  // Compute station response for station 7
  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(kCoordSystem);
  const size_t buffer_size = grid_response->GetStationBufferSize(1);
  std::vector<std::complex<float>> full_response(buffer_size);
  grid_response->Response(BeamMode::kFull, full_response.data(), kTime,
                          kFrequency, station_index, 0);

  // Compute point response
  std::unique_ptr<PointResponse> point_response =
      telescope->GetPointResponse(kTime);
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

BOOST_AUTO_TEST_CASE(under_horizon_point) {
  const Options options;
  std::unique_ptr<Telescope> telescope = LoadAartfaacTelescope(options);

  const size_t station_index = 0;

  // Compute point response
  std::unique_ptr<PointResponse> point_response =
      telescope->GetPointResponse(kTime);
  aocommon::MC2x2F value = aocommon::MC2x2F::Unity();
  point_response->Response(BeamMode::kFull, value.Data(), kCoordSystem.ra,
                           -M_PI /* always under the horizon for Aartfaac */,
                           kFrequency, station_index, 0);

  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_EQUAL(value[i].real(), 0.0);
    BOOST_CHECK_EQUAL(value[i].imag(), 0.0);
  }
}

BOOST_AUTO_TEST_CASE(under_horizon_grid) {
  const Options options;
  std::unique_ptr<Telescope> telescope = LoadAartfaacTelescope(options);
  const size_t width = 16;
  const size_t height = 16;
  const CoordinateSystem under_horizon = {width, height, kRa,     kDec,
                                          0.4,   0.4,    kShiftL, kShiftM};
  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(under_horizon);
  std::vector<std::complex<float>> response(4 * width * height);
  grid_response->Response(everybeam::BeamMode::kElement, response.data(), kTime,
                          kFrequency, 0, 0);

  for (size_t i = 0; i != 4 * width * height; ++i) {
    BOOST_CHECK(std::isfinite(response[i].real()));
    BOOST_CHECK(std::isfinite(response[i].imag()));
  }
  // Check that the top left corner is zero
  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK_EQUAL(response[i].real(), 0.0);
    BOOST_CHECK_EQUAL(response[i].imag(), 0.0);
  }
}

static std::unique_ptr<everybeam::telescope::Telescope> GetTelescopeLofar() {
  Options options;
  options.element_response_model = ElementResponseModel::kLOBES;
  options.coeff_path = LOBES_COEFF_PATH;

  casacore::MeasurementSet ms(LOFAR_LBA_MOCK_MS);
  std::unique_ptr<everybeam::telescope::Telescope> telescope =
      Load(ms, options);

  // AARTFAAC's element response is validated against LOFAR's response. Make
  // sure the needed stations are at the expected offset. When this fails this
  // fix also needs to be applied to the point_response_lobes test.
  const LOFAR& lofar = static_cast<const LOFAR&>(*telescope.get());
  BOOST_REQUIRE_EQUAL(lofar.GetStation(1).GetName(), "CS002LBA");
  BOOST_REQUIRE_EQUAL(lofar.GetStation(2).GetName(), "CS003LBA");
  BOOST_REQUIRE_EQUAL(lofar.GetStation(3).GetName(), "CS004LBA");
  BOOST_REQUIRE_EQUAL(lofar.GetStation(4).GetName(), "CS005LBA");
  BOOST_REQUIRE_EQUAL(lofar.GetStation(5).GetName(), "CS006LBA");
  BOOST_REQUIRE_EQUAL(lofar.GetStation(6).GetName(), "CS007LBA");

  return telescope;
}

static std::unique_ptr<everybeam::telescope::Telescope> GetTelescopeAartfaac() {
  Options options;
  options.element_response_model = ElementResponseModel::kLOBES;
  options.coeff_path = LOBES_COEFF_PATH;
  return LoadAartfaac6Telescope(options);
}

BOOST_AUTO_TEST_CASE(point_response_lobes, *boost::unit_test::tolerance(1e-8)) {
  // Initialize the LOFAR response.
  std::unique_ptr<everybeam::telescope::Telescope> telescope_lofar =
      GetTelescopeLofar();

  std::unique_ptr<everybeam::pointresponse::PointResponse> response_lofar =
      telescope_lofar->GetPointResponse(kTime);

  auto& point_lofar = dynamic_cast<everybeam::pointresponse::PhasedArrayPoint&>(
      *response_lofar);

  // Initialize the AARTFAAC response.
  std::unique_ptr<everybeam::telescope::Telescope> telescope_aartfaac =
      GetTelescopeAartfaac();

  std::unique_ptr<everybeam::pointresponse::PointResponse> response_aartfaac =
      telescope_aartfaac->GetPointResponse(kTime);

  auto& point_aartfaac =
      dynamic_cast<everybeam::pointresponse::PhasedArrayPoint&>(
          *response_aartfaac);

  // Since the coordinate systems of LOFAR and AARTFAAC aren't identical the
  // comparison of the responses needs to be done in the local coordinate
  // system.
  point_lofar.SetUseLocalCoordinateSystem(true);
  point_aartfaac.SetUseLocalCoordinateSystem(true);

  for (int i = 0; i < 288; ++i) {
    // For LOFAR has 6 stations with 48 elements each.
    // For AAFTFAAC has 288 stations with 1 element each.
    // The core stations used on AARTFAAC are LOFAR stations [1, 7].

    const int lofar_station = 1 + i / 48;
    const int lofar_element = i % 48;

    const aocommon::MC2x2 response_lofar = point_lofar.ElementResponse(
        lofar_station, kFrequency, kDirection, lofar_element);
    const aocommon::MC2x2 response_aartfaac =
        point_aartfaac.ElementResponse(i, kFrequency, kDirection, 0);

    BOOST_CHECK_EQUAL(response_lofar[0], response_aartfaac[0]);
    BOOST_CHECK_EQUAL(response_lofar[1], response_aartfaac[1]);
    BOOST_CHECK_EQUAL(response_lofar[2], response_aartfaac[2]);
    BOOST_CHECK_EQUAL(response_lofar[3], response_aartfaac[3]);
  }
}

BOOST_AUTO_TEST_CASE(element_is_full, *boost::unit_test::tolerance(1e-8)) {
  std::unique_ptr<everybeam::telescope::Telescope> telescope =
      GetTelescopeAartfaac();

  std::unique_ptr<everybeam::pointresponse::PointResponse> response =
      telescope->GetPointResponse(kTime);

  std::array<std::complex<float>, 4> buffer_element;
  std::array<std::complex<float>, 4> buffer_full;
  for (int i = 0; i < 288; ++i) {
    response->Response(everybeam::BeamMode::kElement, buffer_element.data(),
                       kRa, kDec, kFrequency, i, 0);
    response->Response(everybeam::BeamMode::kFull, buffer_full.data(), kRa,
                       kDec, kFrequency, i, 0);

    BOOST_CHECK_EQUAL(buffer_element[0], buffer_full[0]);
    BOOST_CHECK_EQUAL(buffer_element[1], buffer_full[1]);
    BOOST_CHECK_EQUAL(buffer_element[2], buffer_full[2]);
    BOOST_CHECK_EQUAL(buffer_element[3], buffer_full[3]);
  }
}

BOOST_AUTO_TEST_CASE(gridded_response_lobes,
                     *boost::unit_test::tolerance(1e-8)) {
  Options options;
  options.element_response_model = ElementResponseModel::kLOBES;
  std::unique_ptr<Telescope> telescope = LoadAartfaac6Telescope(options);

  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(kCoordSystem);
  BOOST_CHECK(dynamic_cast<AartfaacGrid*>(grid_response.get()));

  const size_t buffer_size = grid_response->GetStationBufferSize(1);
  BOOST_CHECK_EQUAL(buffer_size, 4u * kWidth * kHeight);

  std::vector<std::complex<float>> full_response(buffer_size);
  std::vector<std::complex<float>> array_factor(buffer_size);
  std::vector<std::complex<float>> element_response(buffer_size);

  // Reference solution element response
  const size_t buffer_size_all =
      grid_response->GetStationBufferSize(telescope->GetNrStations());
  std::vector<std::complex<float>> element_response_ref(buffer_size_all);
  grid_response->ResponseAllStations(everybeam::BeamMode::kElement,
                                     element_response_ref.data(), kTime,
                                     kFrequency, 0);

  // Reference solution array factor (identity)
  const std::vector<std::complex<float>> array_factor_ref = {
      {1, 0}, {0, 0}, {0, 0}, {1, 0}};

  // Compute station responses all at once
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
        element_response_ref.data() + i * buffer_size,
        element_response_ref.data() + (i + 1) * buffer_size);

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

BOOST_AUTO_TEST_SUITE_END()
