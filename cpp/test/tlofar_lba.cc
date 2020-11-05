#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "../load.h"
#include "../options.h"
#include "../griddedresponse/lofargrid.h"
#include "../elementresponse.h"
#include "../station.h"
#include "../common/types.h"
#include "../../external/npy.hpp"
#include "../telescope/lofar.h"

#include "config.h"
#include <complex>
#include <cmath>

using everybeam::ElementResponseModel;
using everybeam::Load;
using everybeam::matrix22c_t;
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
  // Only checks if we are defaulting to LOBES. Effectively,
  // all the computations will be done as if the Hamaker model was chosen
  // except for station 20 (CS302LBA)
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
  matrix22c_t target_element_response = {{{0}}};
  target_element_response[0][0] = {-0.802669, 0.00378276};
  target_element_response[0][1] = {-0.577012, 0.000892636};
  target_element_response[1][0] = {-0.586008, 0.00549141};
  target_element_response[1][1] = {0.805793, -0.00504886};

  const Station& station =
      static_cast<const Station&>(*(lofartelescope.GetStation(19).get()));
  matrix22c_t element_response =
      station.ComputeElementResponse(time, frequency, direction, true);

  for (size_t i = 0; i != 2; ++i) {
    for (size_t j = 0; j != 2; ++j) {
      BOOST_CHECK(std::abs(element_response[i][j] -
                           target_element_response[i][j]) < 1e-6);
    }
  }

  // Compute station response for station 31 (see also python/test)
  const Station& station31 =
      static_cast<const Station&>(*(lofartelescope.GetStation(31).get()));

  // Channel frequency of channel 4 (3 given zero-based indexing)
  double freq4 = lofartelescope.GetChannelFrequency(3);
  vector3r_t direction_s31 = {0.667806, -0.0770635, 0.740335};
  vector3r_t station0_dir = {0.655743, -0.0670973, 0.751996};
  matrix22c_t station31_response = station31.Response(
      time, freq4, direction_s31, freq4, station0_dir, station0_dir);

  matrix22c_t target_station_response = {{{0}}};
  target_station_response[0][0] = {-0.71383788, 0.00612506};
  target_station_response[0][1] = {-0.4903527, 0.00171652};
  target_station_response[1][0] = {-0.502122, 0.00821683};
  target_station_response[1][1] = {0.7184408, -0.00821723};

  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      BOOST_CHECK(std::abs(station31_response[i][j] -
                           target_station_response[i][j]) < 1e-6);
    }
  }

  // Gridded response
  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(coord_system);
  BOOST_CHECK(nullptr != dynamic_cast<LOFARGrid*>(grid_response.get()));

  // Define buffer and get gridded responses
  std::vector<std::complex<float>> antenna_buffer_single(
      grid_response->GetStationBufferSize(1));

  grid_response->CalculateStation(antenna_buffer_single.data(), time, frequency,
                                  31, 0);

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

BOOST_AUTO_TEST_CASE(test_lobes) {
  Options options;
  // Effectively, all the computations will be done as if the Hamaker model was
  // chosen except for station 20 (CS302LBA)
  options.element_response_model = ElementResponseModel::kLOBES;

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
  grid_response->CalculateStation(antenna_buffer_single.data(), time, frequency,
                                  20, 0);

  // Compare with everybeam at pixel (1, 3). This solution only is a "reference"
  // certainly not a "ground-truth"
  std::vector<std::complex<float>> everybeam_ref_p13 = {
      {0.022321859, 0.032528818},
      {-0.0022747002, 0.096434057},
      {0.037618704, 0.015968619},
      {-0.041608963, 0.021379814}};
  std::size_t offset_13 = (3 + 1 * width) * 4;

  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK_CLOSE(antenna_buffer_single[offset_13 + i],
                      everybeam_ref_p13[i], 1e-4);
  }
  // const long unsigned leshape[] = {(long unsigned int)width, height, 2, 2};
  // npy::SaveArrayAsNumpy("lobes_station_response.npy", false, 4, leshape,
  //                       antenna_buffer_single);
}

BOOST_AUTO_TEST_SUITE_END()
