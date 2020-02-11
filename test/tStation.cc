#include <iostream>

#include "Station.h"

// #include "BeamFormer.h"
#include "LofarMetaDataUtil.h"
#include "MathUtil.h"

using namespace LOFAR::StationResponse;

int main()
{
    std::string name = "CS001LBA";
    vector3r_t position;
//     ElementResponseModel model = ElementResponseModel::Hamaker;

//     Station station(name, position, model);
//
//     double time;
//     double freq;
//     vector3r_t direction = {0.0, 0.0, 1.0};
//     matrix22c_t response;

//     auto antenna0 = Element::Ptr(new Element(station.get_element_response(),0));
//     auto antenna1 = Element::Ptr(new Element(station.get_element_response(),1));
//
//     station.set_antenna(antenna0);
//
//     std::cout << response[0][0] << std::endl;
//
//     response = station.response(time, freq, direction);
//
//     std::cout << response[0][0] << std::endl;
//
//     auto beam_former = BeamFormer::Ptr(new BeamFormer());
//     beam_former->add_antenna(antenna0);
//     beam_former->add_antenna(antenna1);
//     station.set_antenna(beam_former);
//     response = station.response(time, freq, direction);
//     std::cout << response[0][0] << std::endl;
//

//     casacore::MeasurementSet ms("/home/vdtol/data/imagtest96.MS");
    casacore::MeasurementSet ms("/home/vdtol/imagtest_ion1/imagtest.MS");

    double firstTime = casacore::ScalarColumn<double>(ms, "TIME")(0);
    std::cout << firstTime << std::endl;
    double time = firstTime;
    double freq = 55e6;
    matrix22c_t response;
    vector3r_t direction = {0.0, 0.0, 1.0};


//     Station::Ptr station = readStation(ms, 0, ElementResponseModel::OSKARDipole);
    Station::Ptr station = readStation(ms, 0);

    auto freq_beamformer = freq;
    auto station_pointing = direction;
    auto tile_pointing = direction;

    for(int i=0; i<10; i++) {
        auto d = direction;
        d[1] = -0.2 + 0.04*i;
        d = normalize(d);
        response = station->response(time, freq, d, freq_beamformer, station_pointing, tile_pointing);
        std::cout << response[0][0] << " " << response[0][1] << " " << response[1][0] << " " << response[1][1] << " " << std::endl;
    }

    return 0;
}
