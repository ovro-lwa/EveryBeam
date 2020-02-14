#include "config.h"

#include "tElementBeamCommon.h"

using namespace LOFAR::StationResponse;

int main(int argc, char** argv)
{
    ElementResponseModel model(OSKARSphericalWave);
    double frequency = 140e6; // Mhz
    std::string input_filename(TEST_MEASUREMENTSET);
    std::string output_filename("element-beams-oskar-sphericalwave.fits");
    run(model, frequency, input_filename, output_filename);
}