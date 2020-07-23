#include <cmath>
#include <iostream>
#include <cstdlib>

#include <oskarelementresponse.h>

#include "../../external/npy.hpp"
// #include "npy.hpp"  // to save arrays in numpy format

int main(int argc, char** argv){
// int main() {
    everybeam::OSKARElementResponseSphericalWave element_response("oskar.h5");
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
        double x = (2.0*i)/(N-1) - 1.0;
        for(int j=0; j<N; ++j) {
            double y = (2.0*j)/(N-1) - 1.0;
            double theta = asin(sqrt(x*x + y*y));
            double phi = atan2(y,x);
            element_response.Response(0, freq, theta, phi, result_arr[i][j]);
        }
    }

    const long unsigned leshape [] = {(long unsigned int) N, N, 2, 2};
    npy::SaveArrayAsNumpy("response.npy", false, 4, leshape, result);
}
