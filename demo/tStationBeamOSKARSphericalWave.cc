#include "config.h"

#include "tStationBeamCommon.h"

using namespace everybeam;

int main(int argc, char** argv) {
  ElementResponseModel model(kOSKARSphericalWave);
  double frequency = 140e6;  // Mhz
  std::string input_filename(TEST_MEASUREMENTSET);
  std::string output_filename("station-beams-oskar-sphericalwave.fits");
  run(model, frequency, input_filename, output_filename);
}