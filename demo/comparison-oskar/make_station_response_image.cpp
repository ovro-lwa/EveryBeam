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
#include "telescope/oskar.h"

using namespace everybeam;

int main(int argc, char** argv){

    ElementResponseModel response_model = ElementResponseModel::kOSKARSphericalWave;
    Options options;
    options.element_response_model = response_model;

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
    auto element_response = std::make_shared<everybeam::OSKARElementResponseSphericalWave>("oskar.h5");
    station->SetResponse(element_response);

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

            auto result = station->Response(time, freq, direction, freq0, station0, tile0, false);
            result_arr[i][j][0][0] = result[0][0];
            result_arr[i][j][0][1] = result[0][1];
            result_arr[i][j][1][0] = result[1][0];
            result_arr[i][j][1][1] = result[1][1];

        }
    }

    const long unsigned leshape [] = {(long unsigned int) N, (long unsigned int) N, 2, 2};
    npy::SaveArrayAsNumpy("station-response.npy", false, 4, leshape, result);
}
