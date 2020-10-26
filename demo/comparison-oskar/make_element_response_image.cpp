#include <cmath>
#include <iostream>
#include <cstdlib>

#include <oskarelementresponse.h>

#include "../../external/npy.hpp"   // to save arrays in numpy format

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
    for(int i=0; i<N; ++i) {
        double x = (2.0*i)/(N-1) - 1.0;
        for(int j=0; j<N; ++j) {
            double y = (2.0*j)/(N-1) - 1.0;
            double theta = std::asin(sqrt(x*x + y*y));
            double phi = std::atan2(y,x);
            typedef std::complex<double> Mat[2][2]; 
            element_response.Response(0, freq, theta, phi, *reinterpret_cast<Mat*>(&result[4 * (i*N + j)]));
        }
    }

    const long unsigned leshape [] = {(long unsigned int) N, (long unsigned int) N, 2, 2};
    npy::SaveArrayAsNumpy("response.npy", false, 4, leshape, result);
}
