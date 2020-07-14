#include <boost/test/unit_test.hpp>

#include "./../load.h"
#include "./../options.h"
#include "./../gridded_response/lofargrid.h"
#include "./../element_response.h"
#include "../../external/npy.hpp"

#include "config.h"
#include <complex>
#include <cmath>

using namespace everybeam;

BOOST_AUTO_TEST_CASE(load_lofar) {
  ElementResponseModel response_model = ElementResponseModel::Hamaker;
  Options options;
  casacore::MeasurementSet ms(LOFAR_MOCK_MS);

  // Load LOFAR Telescope
  std::unique_ptr<telescope::Telescope> telescope =
      Load(ms, response_model, options);

  // Assert if we indeed have a LOFAR pointer
  BOOST_CHECK(nullptr != dynamic_cast<telescope::LOFAR*>(telescope.get()));

  // Assert if correct number of stations
  std::size_t nstations = 70;
  BOOST_CHECK_EQUAL(telescope->GetNrStations(), nstations);

  // Assert if GetStation(stationd_id) behaves properly
  BOOST_CHECK_EQUAL(telescope->GetStation(0)->name(), "CS001HBA0");

  // Properties extracted from MS
  double time = 4929192878.008341;
  double frequency = 138476562.5;
  std::size_t width(4), height(4);
  double ra(2.15374123), dec(0.8415521), dl(0.5 * M_PI / 180.),
      dm(0.5 * M_PI / 180.), shift_l(0.), shift_m(0.);

  coords::CoordinateSystem coord_system = {.width = width,
                                           .height = height,
                                           .ra = ra,
                                           .dec = dec,
                                           .dl = dl,
                                           .dm = dm,
                                           .phase_centre_dl = shift_l,
                                           .phase_centre_dm = shift_m};
  std::unique_ptr<gridded_response::GriddedResponse> grid_response =
      telescope->GetGriddedResponse(coord_system);
  BOOST_CHECK(nullptr !=
              dynamic_cast<gridded_response::LOFARGrid*>(grid_response.get()));

  // Define buffer and get gridded responses
  std::vector<std::complex<float>> antenna_buffer_single(
      grid_response->GetBufferSize(1));
  grid_response->CalculateStation(antenna_buffer_single.data(), time, frequency,
                                  23);
  BOOST_CHECK_EQUAL(antenna_buffer_single.size(),
                    std::size_t(width * height * 2 * 2));

  // LOFARBeam output at pixel (2,2):
  std::vector<std::complex<float>> lofar_p22 = {{-0.175908, -0.000478397},
                                                {-0.845988, -0.00121503},
                                                {-0.89047, -0.00125383},
                                                {0.108123, -5.36076e-05}};

  // Compare with everybeam
  std::size_t offset_22 = (2 + 2 * height) * 4;
  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK(std::abs(antenna_buffer_single[offset_22 + i] - lofar_p22[i]) <
                1e-4);
  }

  // LOFARBeam output at pixel (1,3):
  std::vector<std::complex<float>> lofar_p13 = {{-0.158755, -0.000749433},
                                                {-0.816165, -0.00272568},
                                                {-0.863389, -0.00283979},
                                                {0.0936919, 0.000110673}};

  // Compare with everybeam
  std::size_t offset_13 = (1 + 3 * height) * 4;
  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK(std::abs(antenna_buffer_single[offset_13 + i] - lofar_p13[i]) <
                1e-4);
  }

  std::vector<std::complex<float>> antenna_buffer_all(
      grid_response->GetBufferSize(telescope->GetNrStations()));
  grid_response->CalculateAllStations(antenna_buffer_all.data(), time,
                                      frequency);
  BOOST_CHECK_EQUAL(
      antenna_buffer_all.size(),
      std::size_t(telescope->GetNrStations() * width * height * 2 * 2));

  // Print to np array
  // const long unsigned leshape[] = {(long unsigned int)width, height, 2, 2};
  // npy::SaveArrayAsNumpy("station_responses.npy", false, 4, leshape,
  //                       antenna_buffer_single);
}