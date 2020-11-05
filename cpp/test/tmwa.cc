#include <boost/test/unit_test.hpp>

#include "../load.h"
#include "../options.h"
#include "../griddedresponse/mwagrid.h"
#include "../../external/npy.hpp"
#include "../telescope/mwa.h"

#include "config.h"
#include <complex>
#include <cmath>

using everybeam::Load;
using everybeam::Options;
using everybeam::coords::CoordinateSystem;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::MWAGrid;
using everybeam::telescope::MWA;
using everybeam::telescope::Telescope;

BOOST_AUTO_TEST_SUITE(tmwa)

BOOST_AUTO_TEST_CASE(load_mwa) {
  Options options;
  options.frequency_interpolation = false;
  options.coeff_path = MWA_COEFF_PATH;

  casacore::MeasurementSet ms(MWA_MOCK_MS);

  std::unique_ptr<Telescope> telescope = Load(ms, options);

  // Assert if we indeed have a MWA pointer
  BOOST_CHECK(nullptr != dynamic_cast<MWA*>(telescope.get()));

  // Assert if correct number of stations
  std::size_t nstations = 128;
  BOOST_CHECK_EQUAL(telescope->GetNrStations(), nstations);

  double time = 4.87541808e+09;
  double frequency = 133794999.99999999;
  std::size_t width(16), height(16);
  double ra(2.18166148), dec(-0.74612826), dl(1. * M_PI / 180.),
      dm(1. * M_PI / 180.), shift_l(0.), shift_m(0.);
  CoordinateSystem coord_system = {.width = width,
                                   .height = height,
                                   .ra = ra,
                                   .dec = dec,
                                   .dl = dl,
                                   .dm = dm,
                                   .phase_centre_dl = shift_l,
                                   .phase_centre_dm = shift_m};

  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(coord_system);
  BOOST_CHECK(nullptr != dynamic_cast<MWAGrid*>(grid_response.get()));

  std::vector<std::complex<float>> antenna_buffer(
      grid_response->GetStationBufferSize(telescope->GetNrStations()));

  try {
    grid_response->CalculateAllStations(antenna_buffer.data(), time, frequency,
                                        0);
  } catch (std::exception& e) {
    throw std::runtime_error(
        std::string(e.what()) +
        "\nPath used for MWA coefficients file: " + MWA_COEFF_PATH);
  } catch (...) {
    throw std::runtime_error(
        std::string("Unknown exception type was thrown, probably by HDF5 "
                    "lib.\nPath used for MWA coefficients file: ") +
        MWA_COEFF_PATH);
  }

  // Check whether pixels are reproduced correctly at certain pixels on 16x16
  // image MWABeam output at pixel (0, 8):
  std::vector<std::complex<float>> mwa_p08 = {{-0.2311999, 0.0294384},
                                              {-0.3220813, 0.03972803},
                                              {-0.303604, 0.0423377},
                                              {0.24611373, -0.03317199}};
  // MWABeam output at pixel (10, 13):
  std::vector<std::complex<float>> mwa_p1013 = {{-0.12771748, 0.02835142},
                                                {-0.0854997, 0.01720171},
                                                {-0.07546211, 0.01246299},
                                                {0.13993847, -0.02099818}};
  // MWABeam output at pixel (15, 15):
  std::vector<std::complex<float>> mwa_p1515 = {{-0.0289045, 0.01751916},
                                                {-0.01501623, 0.0077554},
                                                {-0.01122625, 0.00331686},
                                                {0.02796424, -0.00624412}};

  // Convert pixel to buffer offsets
  std::size_t offset_08 = (0 + 8 * width) * 4;
  std::size_t offset_1013 = (10 + 13 * width) * 4;
  std::size_t offset_1515 = (15 + 15 * width) * 4;

  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK(std::abs(antenna_buffer[offset_08 + i] - mwa_p08[i]) < 1e-6);
    BOOST_CHECK(std::abs(antenna_buffer[offset_1013 + i] - mwa_p1013[i]) <
                1e-6);
    BOOST_CHECK(std::abs(antenna_buffer[offset_1515 + i] - mwa_p1515[i]) <
                1e-6);
  }

  // Print to np array
  const long unsigned leshape[] = {(long unsigned int)width, height, 2, 2};
  npy::SaveArrayAsNumpy("mwa_station_responses.npy", false, 4, leshape,
                        antenna_buffer);
}

BOOST_AUTO_TEST_SUITE_END()
