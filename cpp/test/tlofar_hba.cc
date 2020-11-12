#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "../load.h"
#include "../options.h"
#include "../griddedresponse/lofargrid.h"
#include "../elementresponse.h"
#include "../../external/npy.hpp"
#include "../station.h"
#include "../common/types.h"
#include "../telescope/lofar.h"
#include "../aterms/atermconfig.h"
#include "../aterms/parsetprovider.h"

#include "config.h"
#include <complex>
#include <cmath>
#include <iostream>

using everybeam::ATermSettings;
using everybeam::ElementResponseModel;
using everybeam::Load;
using everybeam::matrix22c_t;
using everybeam::Options;
using everybeam::Station;
using everybeam::vector3r_t;
using everybeam::aterms::ATermConfig;
using everybeam::aterms::ParsetProvider;
using everybeam::coords::CoordinateSystem;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::LOFARGrid;
using everybeam::telescope::LOFAR;
using everybeam::telescope::Telescope;

namespace {
// A very simple override of the ParsetProvider. Just falls back on the
// default or on hard-coded values
struct ParsetATerms : public ParsetProvider {
  virtual std::string GetString(const std::string& key) const final override {
    // Not relevant for EveryBeamATerm
    return "";
  }

  std::string GetStringOr(const std::string& key,
                          const std::string& or_value) const final override {
    // Default response model
    return or_value;
  }

  std::vector<std::string> GetStringList(
      const std::string& key) const final override {
    return std::vector<std::string>{"beam"};
  }

  double GetDoubleOr(const std::string& key,
                     double or_value) const final override {
    // Update interval set to 1200 (s)
    return 1200.0;
  }
  bool GetBool(const std::string& key) const final override {
    // No use
    return false;
  }
  bool GetBoolOr(const std::string& key, bool or_value) const final override {
    return or_value;
  }
};
}  // namespace

struct HBAFixture {
  HBAFixture() : time(4929192878.008341), frequency(138476562.5) {
    options.element_response_model = ElementResponseModel::kHamaker;
    ms = casacore::MeasurementSet{LOFAR_HBA_MOCK_MS};
    telescope = Load(ms, options);
    coord_system.width = 4;
    coord_system.height = 4;
    coord_system.ra = 2.15374123;
    coord_system.dec = 0.8415521;
    coord_system.dl = 0.5 * M_PI / 180.;
    coord_system.dm = 0.5 * M_PI / 180.;
    coord_system.phase_centre_dl = 0.;
    coord_system.phase_centre_dm = 0.;
    grid_response = telescope->GetGriddedResponse(coord_system);
  }
  ~HBAFixture(){};
  Options options;
  std::unique_ptr<Telescope> telescope;
  std::unique_ptr<GriddedResponse> grid_response;
  CoordinateSystem coord_system;

  casacore::MeasurementSet ms;
  double time;
  double frequency;
};

BOOST_FIXTURE_TEST_SUITE(tlofar_hba, HBAFixture)

BOOST_AUTO_TEST_CASE(load_lofar_hba) {
  options.element_response_model = ElementResponseModel::kHamaker;
  // Assert if we indeed have a LOFAR pointer
  BOOST_CHECK(nullptr != dynamic_cast<LOFAR*>(telescope.get()));

  // Assert if correct number of stations
  std::size_t nstations = 70;
  BOOST_CHECK_EQUAL(telescope->GetNrStations(), nstations);

  // Assert if GetStation(stationd_id) behaves properly
  const LOFAR& lofartelescope = static_cast<const LOFAR&>(*telescope.get());
  BOOST_CHECK_EQUAL(lofartelescope.GetStation(0)->GetName(), "CS001HBA0");
}

BOOST_AUTO_TEST_CASE(element_response) {
  const LOFAR& lofartelescope = static_cast<const LOFAR&>(*telescope.get());
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
  matrix22c_t target_element_response = {{{0}}};
  target_element_response[0][0] = {-0.164112, -0.000467162};
  target_element_response[0][1] = {-0.843709, -0.00123631};
  target_element_response[1][0] = {-0.892528, -0.00126278};
  target_element_response[1][1] = {0.0968527, -6.7158e-05};

  const Station& station =
      static_cast<const Station&>(*(lofartelescope.GetStation(11).get()));
  matrix22c_t element_response =
      station.ComputeElementResponse(time, frequency, direction, false);

  // Check whether element_response and target_element_response are "equal"
  for (size_t i = 0; i != 2; ++i) {
    for (size_t j = 0; j != 2; ++j) {
      BOOST_CHECK(std::abs(element_response[i][j] -
                           target_element_response[i][j]) < 1e-6);
    }
  }

  // Compute station response for station 63 (see also python/test)
  const Station& station63 =
      static_cast<const Station&>(*(lofartelescope.GetStation(63).get()));

  vector3r_t direction_s63 = {0.424588, 0.4629957, 0.7780411};
  vector3r_t station0_dir = {0.4083262, 0.5273447, 0.7451022};
  vector3r_t tile0_dir = {0.4083268, 0.5273442, 0.7451022};
  matrix22c_t station63_response = station63.Response(
      time, frequency, direction_s63, frequency, station0_dir, tile0_dir);

  matrix22c_t target_station_response = {{{0}}};
  target_station_response[0][0] = {0.032594235, -0.00023045994};
  target_station_response[0][1] = {0.12204097, -0.00091857865};
  target_station_response[1][0] = {0.13063535, -0.0010039175};
  target_station_response[1][1] = {-0.029348446, 0.00023882818};

  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      BOOST_CHECK(std::abs(station63_response[i][j] -
                           target_station_response[i][j]) < 1e-6);
    }
  }
}

BOOST_AUTO_TEST_CASE(gridded_response) {
  BOOST_CHECK(nullptr != dynamic_cast<LOFARGrid*>(grid_response.get()));

  // Define buffer and get gridded responses
  std::vector<std::complex<float>> antenna_buffer_single(
      grid_response->GetStationBufferSize(1));
  grid_response->CalculateStation(antenna_buffer_single.data(), time, frequency,
                                  23, 0);
  BOOST_CHECK_EQUAL(
      antenna_buffer_single.size(),
      std::size_t(coord_system.width * coord_system.height * 2 * 2));
  // LOFARBeam output at pixel (2,2):
  std::vector<std::complex<float>> lofar_p22 = {{-0.175908, -0.000478397},
                                                {-0.845988, -0.00121503},
                                                {-0.89047, -0.00125383},
                                                {0.108123, -5.36076e-05}};

  // Compare with everybeam
  std::size_t offset_22 = (2 + 2 * coord_system.width) * 4;
  for (std::size_t i = 0; i < 4; ++i) {
    // Tolerance is a percentage, so 1e-2 --> 1e-4
    BOOST_CHECK_CLOSE(antenna_buffer_single[offset_22 + i], lofar_p22[i], 1e-2);
  }

  // LOFARBeam output at pixel (1,3):
  std::vector<std::complex<float>> lofar_p13 = {{-0.158755, -0.000749433},
                                                {-0.816165, -0.00272568},
                                                {-0.863389, -0.00283979},
                                                {0.0936919, 0.000110673}};

  // Compare with everybeam
  std::size_t offset_13 = (1 + 3 * coord_system.width) * 4;
  for (std::size_t i = 0; i < 4; ++i) {
    // Tolerance is a percentage, so 1e-2 --> 1e-4
    BOOST_CHECK_CLOSE(antenna_buffer_single[offset_13 + i], lofar_p13[i], 1e-2);
  }

  // All stations
  std::vector<std::complex<float>> antenna_buffer_all(
      grid_response->GetStationBufferSize(telescope->GetNrStations()));
  grid_response->CalculateAllStations(antenna_buffer_all.data(), time,
                                      frequency, 0);
  BOOST_CHECK_EQUAL(antenna_buffer_all.size(),
                    std::size_t(telescope->GetNrStations() *
                                coord_system.width * coord_system.height * 4));

  // Check consistency of values for station 23
  std::size_t offset_s23 = 23 * coord_system.width * coord_system.height * 4;
  for (std::size_t i = 0; i != antenna_buffer_single.size(); ++i) {
    BOOST_CHECK(std::abs(antenna_buffer_all[offset_s23 + i] -
                         antenna_buffer_single[i]) < 1e-6);
  }

  // Check result via aterm calculation
  // Fake the original time for the aterm calculation, accounting for half the
  // update interval
  double time1 = time - 600;
  ATermSettings aterm_settings;
  std::vector<std::complex<float>> aterm_buffer(antenna_buffer_all.size());
  ParsetATerms parset_aterms;
  ATermConfig aterms(telescope->GetNrStations(), coord_system, aterm_settings);
  aterms.Read(ms, parset_aterms);
  aterms.Calculate(aterm_buffer.data(), time1, frequency, 0, nullptr);

  // Check against antenna_buffer_all
  for (std::size_t i = 0; i != aterm_buffer.size(); ++i) {
    BOOST_CHECK_CLOSE(aterm_buffer[i], antenna_buffer_all[i], 1e-6);
  }

  // Save buffer for later reference
  std::vector<std::complex<float>> aterm_ref = aterm_buffer;

  // Result should not change for time increase <1200s
  aterms.Calculate(aterm_buffer.data(), time1 + 1199, frequency, 0, nullptr);
  for (std::size_t i = 0; i != aterm_buffer.size(); ++i) {
    BOOST_CHECK_CLOSE(aterm_buffer[i], aterm_ref[i], 1e-6);
  }

  // Result should change for time increase >=1200s
  aterms.Calculate(aterm_buffer.data(), time1 + 1201, frequency, 0, nullptr);
  for (std::size_t i = 0; i != aterm_buffer.size(); ++i) {
    BOOST_CHECK(std::abs(aterm_buffer[i] - aterm_ref[i]) > 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(differential_beam) {
  // Test with differential beam, single
  Options options_diff_beam = options;
  options_diff_beam.element_response_model = ElementResponseModel::kHamaker;
  options_diff_beam.use_differential_beam = true;

  // Load (a new) LOFAR Telescope
  std::unique_ptr<Telescope> telescope_diff_beam = Load(ms, options_diff_beam);

  std::unique_ptr<GriddedResponse> grid_response_diff_beam =
      telescope_diff_beam->GetGriddedResponse(coord_system);

  std::vector<std::complex<float>> antenna_buffer_diff_beam(
      grid_response_diff_beam->GetStationBufferSize(1));
  grid_response_diff_beam->CalculateStation(antenna_buffer_diff_beam.data(),
                                            time, frequency, 15, 0);

  std::size_t offset_22 = (2 + 2 * coord_system.width) * 4;
  double norm_jones_mat = 0.;
  for (std::size_t i = 0; i < 4; ++i) {
    norm_jones_mat += std::norm(antenna_buffer_diff_beam[offset_22 + i]);
  }
  BOOST_CHECK(std::abs(norm_jones_mat - 2.) < 1e-6);
}

BOOST_AUTO_TEST_CASE(integrated_beam) {
  // Just check whether CalculateIntegratedResponse does run and reproduces
  // results for One time interval
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
  // stationGains[a1].HermTranspose().Transpose(),
  //                                     stationGains[a2]);
  // -          double w = weights.Value(a1, a2);
  // +          double w = 1.;
  //
  std::size_t width_pb = 40, height_pb = 40;
  // (0.25 * M_PI / 180.) equals 900asec
  CoordinateSystem coord_system_pb = {width_pb,
                                      height_pb,
                                      coord_system.ra,
                                      coord_system.dec,
                                      (0.25 * M_PI / 180.),
                                      (0.25 * M_PI / 180.),
                                      coord_system.phase_centre_dl,
                                      coord_system.phase_centre_dm};

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