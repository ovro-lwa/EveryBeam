#include "config.h"

#include "tStationBeamCommon.h"

using namespace LOFAR::StationResponse;

int main(int argc, char** argv)
{
    ElementResponseModel model(OSKARDipole);
    double frequency = 140e6; // Mhz
    std::string input_filename(TEST_MEASUREMENTSET);
    std::string output_filename("station-beams-oskar-dipole.fits");
    run(model, frequency, input_filename, output_filename);
}