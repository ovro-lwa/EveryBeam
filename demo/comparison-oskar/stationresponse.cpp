#include <cmath>
#include <iostream>
#include <cstdlib>

#include <oskarelementresponse.h>

#include "../../external/npy.hpp"   // to save arrays in numpy format

#include "load.h"
#include "options.h"
#include "config.h"
#include "coords/coordutils.h"
#include "coords/itrfconverter.h"

using namespace everybeam;

int main(int argc, char** argv){

    ElementResponseModel response_model = ElementResponseModel::kOSKARSphericalWave;
    Options options;
    options.element_response_model = response_model;

//     casacore::MeasurementSet ms("/home/vdtol/skalowmini/skalowmini-coef1.MS");
    casacore::MeasurementSet ms("skalowmini-coef.MS");

    casacore::ScalarMeasColumn<casacore::MDirection> referenceDirColumn(
        ms.field(),
        casacore::MSField::columnName(casacore::MSFieldEnums::REFERENCE_DIR));

    auto reference_dir = referenceDirColumn(0);

    std::cout << ms.nrow() << std::endl;;
    auto unique_times_table = ms.sort("TIME", casacore::Sort::Ascending, casacore::Sort::NoDuplicates);
    std::cout << unique_times_table.nrow() << std::endl;

    casacore::ScalarColumn<double> time_column(ms, "TIME");

    real_t time = time_column(0);
    std::cout << "time: " << time << std::endl;

    std::cout << reference_dir << std::endl;

    vector3r_t station0 ;
    vector3r_t tile0;
    coords::ITRFConverter itrf_converter(time);
    coords::SetITRFVector(itrf_converter.ToDirection(reference_dir), station0);
    coords::SetITRFVector(itrf_converter.ToDirection(reference_dir), tile0);

    // Load OSKAR Telescope
    auto telescope = Load(ms, options);
    auto &oskar_telescope = *dynamic_cast<telescope::OSKAR*>(telescope.get());

    Station::Ptr station = oskar_telescope.GetStation(0);

    Antenna::Ptr antenna = station->GetAntenna();

    auto p = antenna->coordinate_system_.axes.p;
    auto q = antenna->coordinate_system_.axes.q;
    auto r = antenna->coordinate_system_.axes.r;


    double freq = 50e6;

    int N;
    if (argc == 1){
        N = 256;
    }
    else{
        N = atoi(argv[1]);
    }

    std::vector<std::complex<double>> result(N*N*2*2);
    typedef std::complex<double>result_arr_t[N][N][2][2];
    result_arr_t &result_arr = * (result_arr_t*) result.data();


    for(int i=0; i<N; ++i) {
        std::cout << i << std::endl;
        double x = (2.0*i)/(N-1) - 1.0;
        for(int j=0; j<N; ++j) {
            double y = (2.0*j)/(N-1) - 1.0;

            double z = sqrt(1.0 - x*x - y*y);

            const vector3r_t direction = {
                x*p[0] + y*q[0] + z*r[0],
                x*p[1] + y*q[1] + z*r[1],
                x*p[2] + y*q[2] + z*r[2]
            };
            real_t freq0 = 50e6;

            auto result = station->Response(time, freq, direction, freq0, station0, tile0);
            result_arr[i][j][0][0] = result[0][0];
            result_arr[i][j][0][1] = result[0][1];
            result_arr[i][j][1][0] = result[1][0];
            result_arr[i][j][1][1] = result[1][1];

//             std::cout << result_arr[i][j][0][0] << ", " << result_arr[i][j][1][0] << std::endl;

//             auto result = station->ArrayFactor(time, freq, direction, freq0, station0, tile0);
//             result_arr[i][j][0][0] = result[0];
//             result_arr[i][j][0][1] = 0.0;
//             result_arr[i][j][1][0] = 0.0;
//             result_arr[i][j][1][1] = result[1];

        }
    }

    const long unsigned leshape [] = {(long unsigned int) N, (long unsigned int) N, 2, 2};
    npy::SaveArrayAsNumpy("station-response.npy", false, 4, leshape, result);
}


/*

#include "../options.h"
#include "../griddedresponse/lofargrid.h"
#include "../elementresponse.h"
#include "../../external/npy.hpp"

#include <complex>
#include <cmath>

using namespace everybeam;

BOOST_AUTO_TEST_CASE(load_lofar) {
  ElementResponseModel response_model = ElementResponseModel::kHamaker;
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
  BOOST_CHECK_EQUAL(telescope->GetStation(0)->GetName(), "CS001HBA0");

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
  std::unique_ptr<griddedresponse::GriddedResponse> grid_response =
      telescope->GetGriddedResponse(coord_system);
  BOOST_CHECK(nullptr !=
              dynamic_cast<griddedresponse::LOFARGrid*>(grid_response.get()));

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

  //   std::vector<std::complex<float>> antenna_buffer_all(
  //       grid_response->GetBufferSize(telescope->GetNrStations()));
  //   grid_response->CalculateAllStations(antenna_buffer_all.data(), time,
  //                                       frequency);
  //   BOOST_CHECK_EQUAL(
  //       antenna_buffer_all.size(),
  //       std::size_t(telescope->GetNrStations() * width * height * 2 * 2));

  // Test with differential beam, single
  Options options_diff_beam;
  options_diff_beam.use_differential_beam = true;

  // Load LOFAR Telescope
  std::unique_ptr<telescope::Telescope> telescope_diff_beam =
      Load(ms, response_model, options_diff_beam);

  std::unique_ptr<griddedresponse::GriddedResponse> grid_response_diff_beam =
      telescope_diff_beam->GetGriddedResponse(coord_system);

  std::vector<std::complex<float>> antenna_buffer_diff_beam(
      grid_response_diff_beam->GetBufferSize(1));
  grid_response_diff_beam->CalculateStation(antenna_buffer_diff_beam.data(),
                                            time, frequency, 15);

  double norm_jones_mat = 0.;
  for (std::size_t i = 0; i < 4; ++i) {
    norm_jones_mat += std::norm(antenna_buffer_diff_beam[offset_22 + i]);
  }
  BOOST_CHECK(std::abs(norm_jones_mat - 2.) < 1e-6);

  // Print to np array
  // const long unsigned leshape[] = {(long unsigned int)width, height, 2, 2};
  // npy::SaveArrayAsNumpy("station_responses.npy", false, 4, leshape,
  //                       antenna_buffer_single);
}

*/
