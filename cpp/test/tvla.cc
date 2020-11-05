#include <boost/test/unit_test.hpp>

#include "../load.h"
#include "../options.h"
#include "../griddedresponse/dishgrid.h"
#include "../elementresponse.h"
#include "../telescope/dish.h"
#include "../../external/npy.hpp"

#include "config.h"
#include <complex>
#include <cmath>

using everybeam::Load;
using everybeam::Options;
using everybeam::coords::CoordinateSystem;
using everybeam::griddedresponse::DishGrid;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::telescope::Dish;
using everybeam::telescope::Telescope;

BOOST_AUTO_TEST_SUITE(tvla)

BOOST_AUTO_TEST_CASE(load_vla) {
  Options options;
  casacore::MeasurementSet ms(VLA_MOCK_MS);

  std::unique_ptr<Telescope> telescope = Load(ms, options);

  // Assert if we indeed have a VLA pointer
  BOOST_CHECK(nullptr != dynamic_cast<Dish*>(telescope.get()));

  // Assert if correct number of stations
  std::size_t nstations = 25;
  BOOST_CHECK_EQUAL(telescope->GetNrStations(), nstations);

  // Time is irrelevant for dish telescope
  double time = 0.5 * (4.90683119e+09 + 4.90684196e+09);
  double frequency = 0.5e+09;  // 1991000000.0;
  std::size_t width(16), height(16);
  double ra(2.62880729), dec(0.02831797), dl(0.125 * M_PI / 180.),
      dm(0.125 * M_PI / 180.), shift_l(0.), shift_m(0.);

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
  BOOST_CHECK(nullptr != dynamic_cast<DishGrid*>(grid_response.get()));

  std::vector<std::complex<float>> antenna_buffer(
      grid_response->GetStationBufferSize(telescope->GetNrStations()));
  grid_response->CalculateAllStations(antenna_buffer.data(), time, frequency,
                                      0);

  // Check that XX/YY are equal
  // Loop over pixels to check that off-diagonals are 0
  // and diagonal entries are equal
  for (std::size_t pxl = 0; pxl < width * height; ++pxl) {
    BOOST_CHECK(std::abs(antenna_buffer[pxl * 4 + 1]) < 1e-8);
    BOOST_CHECK(std::abs(antenna_buffer[pxl * 4 + 2]) < 1e-8);
    BOOST_CHECK(
        std::abs(antenna_buffer[pxl * 4] - antenna_buffer[pxl * 4 + 3]) < 1e-8);
  }

  // Check whether pixels are reproduced correctly at certain pixels on 16x16
  // image VLABeam output at pixel (0, 0):
  std::vector<std::complex<float>> vla_p00 = {
      {0.4687362, 0.}, {0, 0}, {0, 0}, {0.4687362, 0.}};
  // VLABeam output at pixel (2, 3):
  std::vector<std::complex<float>> vla_p23 = {
      {0.65400606, 0.}, {0, 0}, {0, 0}, {0.65400606, 0.}};
  // VLABeam output at pixel (10, 12):
  std::vector<std::complex<float>> vla_p1012 = {
      {0.8672509, 0.}, {0, 0}, {0, 0}, {0.8672509, 0.}};

  // Convert pixel to buffer offsets
  std::size_t offset_00 = (0 + 0 * width) * 4;
  std::size_t offset_23 = (2 + 3 * width) * 4;
  std::size_t offset_1012 = (10 + 12 * width) * 4;

  // Check if results are reproduced
  BOOST_CHECK_EQUAL_COLLECTIONS(antenna_buffer.begin() + offset_00,
                                antenna_buffer.begin() + offset_00 + 4,
                                vla_p00.begin(), vla_p00.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(antenna_buffer.begin() + offset_23,
                                antenna_buffer.begin() + offset_23 + 4,
                                vla_p23.begin(), vla_p23.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(antenna_buffer.begin() + offset_1012,
                                antenna_buffer.begin() + offset_1012 + 4,
                                vla_p1012.begin(), vla_p1012.end());

  // Print to np array, note: this spits out the transposed grid
  const long unsigned leshape[] = {(long unsigned int)width, height, 2, 2};
  npy::SaveArrayAsNumpy("vla_station_responses.npy", false, 4, leshape,
                        antenna_buffer);
}

BOOST_AUTO_TEST_SUITE_END()
