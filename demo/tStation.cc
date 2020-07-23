#include <iostream>

#include "./../cpp/station.h"
#include "./../cpp/lofarreadutils.h"
#include "./../cpp/common/mathutils.h"

#include "config.h"

using namespace everybeam;

int main() {
  std::string name = "CS001LBA";
  vector3r_t position;
  //     ElementResponseModel model = ElementResponseModel::kHamaker;

  //     Station station(name, position, model);
  //
  //     double time;
  //     double freq;
  //     vector3r_t direction = {0.0, 0.0, 1.0};
  //     matrix22c_t response;

  //     auto antenna0 = Element::Ptr(new
  //     Element(station.GetElementResponse(),0)); auto antenna1 =
  //     Element::Ptr(new Element(station.GetElementResponse(),1));
  //
  //     station.SetAntenna(antenna0);
  //
  //     std::cout << response[0][0] << std::endl;
  //
  //     Response = station.Response(time, freq, direction);
  //
  //     std::cout << response[0][0] << std::endl;
  //
  //     auto beam_former = BeamFormer::Ptr(new BeamFormer());
  //     beam_former->AddAntenna(antenna0);
  //     beam_former->AddAntenna(antenna1);
  //     station.SetAntenna(beam_former);
  //     Response = station.Response(time, freq, direction);
  //     std::cout << response[0][0] << std::endl;
  //

  casacore::MeasurementSet ms(TEST_MEASUREMENTSET);

  double firstTime = casacore::ScalarColumn<double>(ms, "TIME")(0);
  std::cout << firstTime << std::endl;
  double time = firstTime;
  double freq = 55e6;
  matrix22c_t response;
  vector3r_t direction = {0.0, 0.0, 1.0};

  //     Station::Ptr station = ReadLofarStation(ms, 0,
  //     ElementResponseModel::kOSKARDipole);
  Station::Ptr station = ReadLofarStation(ms, 0);

  auto freq_beamformer = freq;
  auto station_pointing = direction;
  auto tile_pointing = direction;

  for (int i = 0; i < 10; i++) {
    auto d = direction;
    d[1] = -0.2 + 0.04 * i;
    d = normalize(d);
    response = station->Response(time, freq, d, freq_beamformer,
                                 station_pointing, tile_pointing);
    std::cout << response[0][0] << " " << response[0][1] << " "
              << response[1][0] << " " << response[1][1] << " " << std::endl;
  }

  return 0;
}
