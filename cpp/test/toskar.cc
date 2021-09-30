// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "../load.h"
#include "../options.h"
#include "../griddedresponse/oskargrid.h"
#include "../pointresponse/oskarpoint.h"
#include "../elementresponse.h"
#include "../../external/npy.hpp"
#include "../telescope/oskar.h"

#include "config.h"
#include <complex>
#include <cmath>
#include <iostream>

#include "../station.h"

#include "../oskar/oskarelementresponse.h"

using everybeam::ElementResponseModel;
using everybeam::Load;
using everybeam::Options;
using everybeam::Station;
using everybeam::coords::CoordinateSystem;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::OSKARGrid;
using everybeam::pointresponse::OSKARPoint;
using everybeam::pointresponse::PointResponse;
using everybeam::telescope::OSKAR;
using everybeam::telescope::Telescope;

BOOST_AUTO_TEST_SUITE(toskar)

BOOST_AUTO_TEST_CASE(load_oskar) {
  Options options;
  options.element_response_model = ElementResponseModel::kOSKARSphericalWave;

  casacore::MeasurementSet ms(OSKAR_MOCK_MS);

  // Load LOFAR Telescope
  std::unique_ptr<Telescope> telescope = Load(ms, options);

  // Assert if we indeed have an OSKAR pointer
  BOOST_CHECK(nullptr != dynamic_cast<OSKAR*>(telescope.get()));

  // Assert if correct number of stations
  std::size_t nstations = 30;
  BOOST_CHECK_EQUAL(telescope->GetNrStations(), nstations);

  // Assert if GetStation(stationd_id) behaves properly
  const OSKAR& oskartelescope = static_cast<const OSKAR&>(*telescope.get());
  BOOST_CHECK_EQUAL(oskartelescope.GetStation(0)->GetName(), "s0000");

  // Properties extracted from MS
  double time = 4.45353e+09;
  double frequency = 5.0e+07;
  double ra(0.349066), dec(-0.523599);

  // Image properties
  std::size_t width(16), height(16);
  double dl(0.5 * M_PI / 180.), dm(0.5 * M_PI / 180.), shift_l(0.), shift_m(0.);

  CoordinateSystem coord_system = {.width = width,
                                   .height = height,
                                   .ra = ra,
                                   .dec = dec,
                                   .dl = dl,
                                   .dm = dm,
                                   .phase_centre_dl = shift_l,
                                   .phase_centre_dm = shift_m};
  // Get GriddedResponse pointer
  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(coord_system);
  BOOST_CHECK(nullptr != dynamic_cast<OSKARGrid*>(grid_response.get()));

  // Get PointResponse pointer
  std::unique_ptr<PointResponse> point_response =
      telescope->GetPointResponse(time);
  BOOST_CHECK(nullptr != dynamic_cast<OSKARPoint*>(point_response.get()));

  // Define buffer and get gridded response (all stations)
  std::vector<std::complex<float>> antenna_buffer(
      grid_response->GetStationBufferSize(telescope->GetNrStations()));

  grid_response->FullResponseAllStations(antenna_buffer.data(), time, frequency,
                                         0);

  BOOST_CHECK_EQUAL(
      antenna_buffer.size(),
      std::size_t(telescope->GetNrStations() * width * height * 2 * 2));

  // Offset for pixels (8, 1), (8, 8), (5, 13) to check
  std::size_t offset_p81 = (1 + 8 * width) * 4,
              offset_p88 = (8 + 8 * width) * 4,
              offset_p513 = (13 + 5 * width) * 4;

  // Pixel (1, 8), values to be reproduced are
  std::vector<std::complex<float>> oskar_p81 = {{-0.000770827, 0.00734416},
                                                {0.000108162, 0.000871187},
                                                {-0.000259571, 0.00107524},
                                                {0.000648075, -0.00719359}};
  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK(std::abs(antenna_buffer[offset_p81 + i] - oskar_p81[i]) < 1e-6);
  }

  // Pixel (8, 8), values to be reproduced are
  std::vector<std::complex<float>> oskar_p88 = {{-0.00115441, 0.00881435},
                                                {6.07546e-05, 0.0013592},
                                                {-0.000379462, 0.00158725},
                                                {0.000982019, -0.00856565}};
  // Compute results via PointResponse (station 0)
  std::complex<float> point_buffer_single_station[4];
  point_response->FullResponse(point_buffer_single_station, coord_system.ra,
                               coord_system.dec, frequency, 0, 0);

  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK(std::abs(antenna_buffer[offset_p88 + i] - oskar_p88[i]) < 1e-6);
    BOOST_CHECK(std::abs(point_buffer_single_station[i] -
                         antenna_buffer[offset_p88 + i]) < 1e-6);
  }

  // Pixel (5, 13), values to be reproduced are
  std::vector<std::complex<float>> oskar_p513 = {{-0.00114459, 0.00757855},
                                                 {-4.72324e-06, 0.00138143},
                                                 {-0.000382852, 0.00158802},
                                                 {0.000975401, -0.00731452}};
  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK(std::abs(antenna_buffer[offset_p513 + i] - oskar_p513[i]) <
                1e-6);
  }

  // Define buffer and get gridded response (all stations)
  std::vector<std::complex<float>> antenna_buffer_single(
      grid_response->GetStationBufferSize(1));

  // Get FullResponse response from station 5
  grid_response->FullResponse(antenna_buffer_single.data(), time, frequency, 5,
                              0);
  // Check that results in antenna_buffer_single matches the relevant part in
  // antenna_buffer
  std::size_t offset_s5 = 5 * width * height * 4;
  for (std::size_t i = 0; i != antenna_buffer_single.size(); ++i) {
    BOOST_CHECK(std::abs(antenna_buffer[offset_s5 + i] -
                         antenna_buffer_single[i]) < 1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()
