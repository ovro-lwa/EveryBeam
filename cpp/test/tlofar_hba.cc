#include <boost/test/unit_test.hpp>

#include "../load.h"
#include "../options.h"
#include "../griddedresponse/lofargrid.h"
#include "../elementresponse.h"
#include "../../external/npy.hpp"
#include "../station.h"
#include "../common/types.h"

#include "config.h"
#include <complex>
#include <cmath>
#include <iostream>

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

BOOST_AUTO_TEST_SUITE(tlofar_hba)

BOOST_AUTO_TEST_CASE(load_lofar) {
  Options options;
  options.element_response_model = ElementResponseModel::kHamaker;

  casacore::MeasurementSet ms(LOFAR_HBA_MOCK_MS);

  // Load LOFAR Telescope
  std::unique_ptr<Telescope> telescope = Load(ms, options);

  // Assert if we indeed have a LOFAR pointer
  BOOST_CHECK(nullptr != dynamic_cast<LOFAR*>(telescope.get()));

  // Assert if correct number of stations
  std::size_t nstations = 70;
  BOOST_CHECK_EQUAL(telescope->GetNrStations(), nstations);

  // Assert if GetStation(stationd_id) behaves properly
  const LOFAR& lofartelescope = static_cast<const LOFAR&>(*telescope.get());
  BOOST_CHECK_EQUAL(lofartelescope.GetStation(0)->GetName(), "CS001HBA0");

  // Properties extracted from MS
  double time = 4929192878.008341;
  double frequency = 138476562.5;
  std::size_t width(4), height(4);
  double ra(2.15374123), dec(0.8415521), dl(0.5 * M_PI / 180.),
      dm(0.5 * M_PI / 180.), shift_l(0.), shift_m(0.);

  CoordinateSystem coord_system = {.width = width,
                                   .height = height,
                                   .ra = ra,
                                   .dec = dec,
                                   .dl = dl,
                                   .dm = dm,
                                   .phase_centre_dl = shift_l,
                                   .phase_centre_dm = shift_m};
  // Compute and check the Station::ComputeElementResponse for
  // a "randomly selected"  station (station 11)
  // NOTE: this is a regression test in the sense that we only check whether
  // results are consistently reproduced. Reference solution obtained at
  // commit sha 70a286e7dace4616417b0e973a624477f15c9ce3
  //
  // Direction corresponds to one of the itrf directions of the (16) pixels
  // target_element_response is the element response corresponding to this
  // direction
  vector3r_t direction = {0.397408, 0.527527, 0.750855};
  matrix22c_t target_element_response = {0};
  target_element_response[0][0] = {-0.164112, -0.000467162};
  target_element_response[0][1] = {-0.843709, -0.00123631};
  target_element_response[1][0] = {-0.892528, -0.00126278};
  target_element_response[1][1] = {0.0968527, -6.7158e-05};

  const Station& station =
      static_cast<const Station&>(*(lofartelescope.GetStation(11).get()));
  matrix22c_t element_response =
      station.ComputeElementResponse(time, frequency, direction, true);

  // Check whether element_response and target_element_response are "equal"
  for (size_t i = 0; i != 2; ++i) {
    for (size_t j = 0; j != 2; ++j) {
      BOOST_CHECK(std::abs(element_response[i][j] -
                           target_element_response[i][j]) < 1e-6);
    }
  }

  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(coord_system);
  BOOST_CHECK(nullptr != dynamic_cast<LOFARGrid*>(grid_response.get()));

  // Define buffer and get gridded responses
  std::vector<std::complex<float>> antenna_buffer_single(
      grid_response->GetStationBufferSize(1));
  grid_response->CalculateStation(antenna_buffer_single.data(), time, frequency,
                                  23, 0);
  BOOST_CHECK_EQUAL(antenna_buffer_single.size(),
                    std::size_t(width * height * 2 * 2));

  // LOFARBeam output at pixel (2,2):
  std::vector<std::complex<float>> lofar_p22 = {{-0.175908, -0.000478397},
                                                {-0.845988, -0.00121503},
                                                {-0.89047, -0.00125383},
                                                {0.108123, -5.36076e-05}};

  // Compare with everybeam
  std::size_t offset_22 = (2 + 2 * width) * 4;
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
  std::size_t offset_13 = (1 + 3 * width) * 4;
  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK(std::abs(antenna_buffer_single[offset_13 + i] - lofar_p13[i]) <
                1e-4);
  }

  // All stations
  std::vector<std::complex<float>> antenna_buffer_all(
      grid_response->GetStationBufferSize(telescope->GetNrStations()));
  grid_response->CalculateAllStations(antenna_buffer_all.data(), time,
                                      frequency, 0);
  BOOST_CHECK_EQUAL(
      antenna_buffer_all.size(),
      std::size_t(telescope->GetNrStations() * width * height * 2 * 2));

  // Test with differential beam, single
  Options options_diff_beam;
  options_diff_beam.element_response_model = ElementResponseModel::kHamaker;
  options_diff_beam.use_differential_beam = true;

  // Load LOFAR Telescope
  std::unique_ptr<Telescope> telescope_diff_beam = Load(ms, options_diff_beam);

  std::unique_ptr<GriddedResponse> grid_response_diff_beam =
      telescope_diff_beam->GetGriddedResponse(coord_system);

  std::vector<std::complex<float>> antenna_buffer_diff_beam(
      grid_response_diff_beam->GetStationBufferSize(1));
  grid_response_diff_beam->CalculateStation(antenna_buffer_diff_beam.data(),
                                            time, frequency, 15, 0);

  double norm_jones_mat = 0.;
  for (std::size_t i = 0; i < 4; ++i) {
    norm_jones_mat += std::norm(antenna_buffer_diff_beam[offset_22 + i]);
  }
  BOOST_CHECK(std::abs(norm_jones_mat - 2.) < 1e-6);

  // Primary beam tests
  //
  // Just check whether CalculateIntegratedResponse does run and reproduces
  // results One time interval
  std::vector<double> antenna_buffer_integrated(
      grid_response->GetIntegratedBufferSize());
  std::vector<double> baseline_weights(
      telescope->GetNrStations() * (telescope->GetNrStations() + 1) / 2, 1.);
  grid_response->CalculateIntegratedResponse(antenna_buffer_integrated.data(),
                                             time, frequency, 0, 2,
                                             baseline_weights);

  // Just check whether some (rather arbitrary) numbers are reproduced
  BOOST_CHECK(std::abs(antenna_buffer_integrated[10] - 0.0262708) < 1e-6);
  BOOST_CHECK(std::abs(antenna_buffer_integrated[20] - 0.127972) < 1e-6);
  BOOST_CHECK(std::abs(antenna_buffer_integrated.back() - 0.00847742) < 1e-6);

  // Two time intervals, should give same output as single time interval
  std::fill(antenna_buffer_integrated.begin(), antenna_buffer_integrated.end(),
            0);
  std::vector<double> tarray = {time, time};
  baseline_weights.resize(baseline_weights.size() * tarray.size());
  std::fill(baseline_weights.begin(), baseline_weights.end(), 1.);
  grid_response->CalculateIntegratedResponse(antenna_buffer_integrated.data(),
                                             tarray, frequency, 0, 2,
                                             baseline_weights);

  BOOST_CHECK(std::abs(antenna_buffer_integrated[10] - 0.0262708) < 1e-6);
  BOOST_CHECK(std::abs(antenna_buffer_integrated[20] - 0.127972) < 1e-6);
  BOOST_CHECK(std::abs(antenna_buffer_integrated.back() - 0.00847742) < 1e-6);

  // Primary beam response on 40 x 40 grid.
  // Validated results were obtained with the following wsclean command
  //
  // wsclean -size 40 40  -scale 900asec -apply-primary-beam  LOFAR_MOCK.ms
  //
  // where LOFAR_MOCK.ms the MS available from
  // https://www.astron.nl/citt/EveryBeam/L258627-one-timestep.tar.bz2
  //
  // PLEASE NOTE: for the sake of testing, the baseline weights were set to 1
  // in wsclean::lbeamimagemaker, i.e.
  //
  // --- a/lofar/lbeamimagemaker.cpp
  // +++ b/lofar/lbeamimagemaker.cpp
  // @@ -380,7 +380,7 @@ void LBeamImageMaker::makeBeamSnapshot(
  //                                 MC4x4::KroneckerProduct(
  //                                     stationGains[a1].HermTranspose().Transpose(),
  //                                     stationGains[a2]);
  // -          double w = weights.Value(a1, a2);
  // +          double w = 1.;
  //
  std::size_t width_pb = 40, height_pb = 40;
  CoordinateSystem coord_system_pb = {.width = width_pb,
                                      .height = height_pb,
                                      .ra = ra,
                                      .dec = dec,
                                      // 900asec
                                      .dl = (0.25 * M_PI / 180.),
                                      .dm = (0.25 * M_PI / 180.),
                                      .phase_centre_dl = shift_l,
                                      .phase_centre_dm = shift_m};
  std::unique_ptr<GriddedResponse> grid_response_pb =
      telescope->GetGriddedResponse(coord_system_pb);

  std::vector<double> antenna_buffer_pb(
      grid_response_pb->GetIntegratedBufferSize());
  std::vector<double> baseline_weights_pb(
      telescope->GetNrStations() * (telescope->GetNrStations() + 1) / 2, 1.);
  grid_response_pb->CalculateIntegratedResponse(
      antenna_buffer_pb.data(), time, frequency, 0, 8, baseline_weights_pb);
  // Check diagonal and off-diagonal term in component 0 and 5 of HMC4x4
  // representation of Mueller matrix
  std::size_t offset_01616 = 16 * width_pb + 16,
              offset_02310 = 23 * width_pb + 10,
              offset_52020 = 5 * width_pb * height_pb + 20 * width_pb + 20,
              offset_51825 = 5 * width_pb * height_pb + 18 * width_pb + 25;

  BOOST_CHECK(std::abs(antenna_buffer_pb[offset_01616] - 0.020324793) < 1e-6);
  BOOST_CHECK(std::abs(antenna_buffer_pb[offset_02310] - 0.0059926948) < 1e-6);
  BOOST_CHECK(std::abs(antenna_buffer_pb[offset_52020] - 0.00018088287) < 1e-6);
  BOOST_CHECK(std::abs(antenna_buffer_pb[offset_51825] - 0.00013052078) < 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
