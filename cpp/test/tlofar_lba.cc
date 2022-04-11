// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "config.h"
#include "../load.h"
#include "../beammode.h"
#include "../options.h"
#include "../griddedresponse/lofargrid.h"
#include "../elementresponse.h"
#include "../station.h"
#include "../common/types.h"
#include "../telescope/lofar.h"

#include <complex>
#include <cmath>

#include <aocommon/matrix2x2.h>

using everybeam::ElementResponseModel;
using everybeam::Load;
using everybeam::Options;
using everybeam::Station;
using everybeam::vector3r_t;
using everybeam::coords::CoordinateSystem;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::LOFARGrid;
using everybeam::telescope::LOFAR;
using everybeam::telescope::Telescope;

BOOST_AUTO_TEST_SUITE(tlofar_lba)

// Properties extracted from MS
double time = 4.92183348e+09;
double frequency = 57812500.;
double ra(-1.44194878), dec(0.85078091);

// Properties of grid
std::size_t width(4), height(4);
double dl(0.5 * M_PI / 180.), dm(0.5 * M_PI / 180.), shift_l(0.), shift_m(0.);

CoordinateSystem coord_system = {width, height, ra,      dec,
                                 dl,    dm,     shift_l, shift_m};

BOOST_AUTO_TEST_CASE(test_hamaker) {
  Options options;
  options.element_response_model = ElementResponseModel::kHamaker;

  // Load LOFAR Telescope
  std::unique_ptr<Telescope> telescope = Load(LOFAR_LBA_MOCK_MS, options);

  // Assert if we indeed have a LOFAR pointer
  BOOST_CHECK(nullptr != dynamic_cast<LOFAR*>(telescope.get()));
  // Check if
  BOOST_CHECK_EQUAL(telescope->GetNrStations(), size_t{37});

  const LOFAR& lofartelescope = static_cast<const LOFAR&>(*telescope.get());
  BOOST_CHECK_EQUAL(lofartelescope.GetStation(0)->GetName(), "CS001LBA");

  // Reference solution obtained with commit sha
  // 70a286e7dace4616417b0e973a624477f15c9ce3
  //
  // Compute element response for station 19
  // Direction corresponds to the ITRF direction of one of the pixels
  vector3r_t direction = {0.663096, -0.0590573, 0.746199};
  aocommon::MC2x2 target_element_response(
      {-0.802669, 0.00378276}, {-0.577012, 0.000892636},
      {-0.586008, 0.00549141}, {0.805793, -0.00504886});

  const Station& station =
      static_cast<const Station&>(*(lofartelescope.GetStation(19).get()));
  aocommon::MC2x2 element_response =
      station.ComputeElementResponse(time, frequency, direction, false, true);

  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK(std::abs(element_response[i] - target_element_response[i]) <
                1e-6);
  }

  // Compute station response for station 31 (see also python/test)
  const Station& station31 =
      static_cast<const Station&>(*(lofartelescope.GetStation(31).get()));

  // Channel frequency of channel 4 (3 given zero-based indexing)
  double freq4 = lofartelescope.GetChannelFrequency(3);
  vector3r_t direction_s31 = {0.667806, -0.0770635, 0.740335};
  vector3r_t station0_dir = {0.655743, -0.0670973, 0.751996};
  aocommon::MC2x2 station31_response = station31.Response(
      time, freq4, direction_s31, freq4, station0_dir, station0_dir);

  aocommon::MC2x2 target_station_response(
      {-0.71383788, 0.00612506}, {-0.4903527, 0.00171652},
      {-0.502122, 0.00821683}, {0.7184408, -0.00821723});

  for (size_t i = 0; i != 4; ++i) {
    BOOST_CHECK(std::abs(station31_response[i] - target_station_response[i]) <
                1e-6);
  }

  // Gridded response
  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(coord_system);
  BOOST_CHECK(nullptr != dynamic_cast<LOFARGrid*>(grid_response.get()));

  // Define buffer and get gridded responses
  std::vector<std::complex<float>> antenna_buffer_single(
      grid_response->GetStationBufferSize(1));

  grid_response->Response(everybeam::BeamMode::kFull,
                          antenna_buffer_single.data(), time, frequency, 31, 0);

  // Compare with everybeam at pixel (1, 3), reference solution obtained with
  // everybeam at commit sha 70a286e7dace4616417b0e973a624477f15c9ce3
  std::vector<std::complex<float>> everybeam_ref_p13 = {
      {-0.77199, 0.00214581},
      {-0.529293, -0.00127649},
      {-0.541472, 0.0053149},
      {0.775696, -0.00369905}};

  std::size_t offset_13 = (3 + 1 * width) * 4;
  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK(std::abs(antenna_buffer_single[offset_13 + i] -
                         everybeam_ref_p13[i]) < 1e-6);
  }
}

// DOWNLOAD_LOBES required since we need to make sure
// coefficient file for CS302LBA is present in build dir.
#ifdef DOWNLOAD_LOBES
BOOST_AUTO_TEST_CASE(test_lobes) {
  Options options;
  // Effectively, all the computations will be done as if the Hamaker model was
  // chosen except for station 20 (CS302LBA)
  options.element_response_model = ElementResponseModel::kLOBES;
  options.coeff_path = LOBES_COEFF_PATH;

  casacore::MeasurementSet ms(LOFAR_LBA_MOCK_MS);

  // Load LOFAR Telescope
  std::unique_ptr<Telescope> telescope = Load(ms, options);

  // Extract Station 20, should be station CS302LBA
  const LOFAR& lofartelescope = static_cast<const LOFAR&>(*telescope.get());
  BOOST_CHECK_EQUAL(lofartelescope.GetStation(20)->GetName(), "CS302LBA");

  // Gridded response
  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(coord_system);

  // Define buffer and get gridded responses
  std::vector<std::complex<float>> antenna_buffer_single(
      grid_response->GetStationBufferSize(1));

  // Get the gridded response for station 20 (of course!)
  grid_response->Response(everybeam::BeamMode::kFull,
                          antenna_buffer_single.data(), time, frequency, 20, 0);

  // Compare with everybeam at pixel (1, 3). This solution only is a "reference"
  // certainly not a "ground-truth"
  std::vector<std::complex<float>> everybeam_ref_p13 = {
      {-2.1952772, -2.7803752},
      {-1.4990979, -1.9227437},
      {0.0086051673, 0.034685627},
      {1.5801408, 2.0045171}};

  std::size_t offset_13 = (3 + 1 * width) * 4;
  for (std::size_t i = 0; i < 4; ++i) {
    // relative tolerance is 2e-1%, so 2e-3
    BOOST_CHECK_CLOSE(antenna_buffer_single[offset_13 + i],
                      everybeam_ref_p13[i], 2e-1);
  }
}
#endif  // DOWNLOAD_LOBES

BOOST_AUTO_TEST_SUITE_END()
