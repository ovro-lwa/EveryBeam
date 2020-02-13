#include <iostream>
#include <complex>
#include <vector>

#include "config.h"
#include "beam-helper.h"

void calculateStationBeams(
	std::vector<LOFAR::StationResponse::Station::Ptr>& stations,
    std::vector<vector3r_t>& itrfsDirections,
    vector3r_t stationDirection,
    vector3r_t tileDirection,
    bool useDifferentialBeam,
    vector3r_t diffBeamCentre,
    unsigned int subgrid_size,
    std::vector<std::complex<float>>& buffer,
    double time, double frequency)
{
    typedef std::complex<float> Data[stations.size()][subgrid_size][subgrid_size][4];
    Data* data_ptr = (Data *) buffer.data();

    #pragma omp parallel for
	for (size_t s = 0; s < stations.size(); s++) {
        matrix22c_t inverseCentralGain =
            stations[s]->response(time, frequency, diffBeamCentre, frequency, stationDirection, tileDirection);

        for (unsigned y = 0; y < subgrid_size; y++) {
            for (unsigned x = 0; x < subgrid_size; x++) {
                auto direction = itrfsDirections[y * subgrid_size + x];
                auto freq_beamformer = frequency;
                matrix22c_t gainMatrix =
                    stations[s]->response(
                    time, frequency, direction, freq_beamformer,
                    stationDirection, tileDirection);

                std::complex<float>* antBufferPtr = (*data_ptr)[s][y][x];

                matrix22c_t stationGain =
                    useDifferentialBeam ?
                    inverseCentralGain * gainMatrix :
                    gainMatrix;

                antBufferPtr[0] = stationGain[0][0];
                antBufferPtr[1] = stationGain[0][1];
                antBufferPtr[2] = stationGain[1][0];
                antBufferPtr[3] = stationGain[1][1];
            }
        }
    }
}

int main(int argc, char** argv) {
    // Open measurement set
    std::string ms_filename = std::string(TEST_MEASUREMENTSET);
    std::cout << ">> Opening measurementset: " << ms_filename << std::endl;
    casacore::MeasurementSet ms(ms_filename);

    // Read centre frequency
    double centreFrequency = 132e6; // Mhz
    std::clog << "Centre frequency: " << centreFrequency * 1e-6 << " Mhz" << std::endl;

    // Read number of stations
    size_t nr_stations = ms.antenna().nrow();
    std::clog << "Number of stations: " << nr_stations << std::endl;

    // Read number of timesteps
    size_t nr_timesteps = ms.nrow();
    std::clog << "Number of timesteps: " << nr_timesteps << std::endl;

    // Read field id
    casacore::ROScalarColumn<int> fieldIdColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
    size_t fieldId = fieldIdColumn.getColumn()[0];
    std::clog << "Field ID: " << fieldId << std::endl;

    // Read phase centre info
    double phaseCentreRA, phaseCentreDec;
    GetPhaseCentreInfo(ms, fieldId, phaseCentreRA, phaseCentreDec);
    std::clog << "RA:  " << phaseCentreRA << std::endl;
    std::clog << "DEC: " << phaseCentreDec << std::endl;

    // Read observation time
    casacore::ScalarColumn<double> timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
    double currentTime = timeColumn(nr_timesteps / 2);

    // Set element responde model
    LOFAR::StationResponse::ElementResponseModel elementResponseModel(LOFAR::StationResponse::Hamaker);

    // Read stations
	std::vector<LOFAR::StationResponse::Station::Ptr> stations;
    stations.resize(nr_stations);
	readStations(ms, stations.begin(), elementResponseModel);

    // Imaging parameters
    float image_size = 0.5; // in radians
    size_t subgrid_size = 32; // in pixels

    // Evaluate beam directions in ITRF coordinates
    std::cout << ">>> Computing directions to evaluate beam" << std::endl;
    std::vector<vector3r_t> itrfDirections(subgrid_size * subgrid_size);
    GetITRFDirections(itrfDirections.data(), subgrid_size, image_size, currentTime, phaseCentreRA, phaseCentreDec);

    // Set station beam direction to centre of field
    vector3r_t stationDirection = itrfDirections[(subgrid_size/2) * subgrid_size + (subgrid_size/2)];

    // Set tile beam direction equal to station direction
    vector3r_t tileDirection = stationDirection;

    // Set differential beam direction equal to station direction
    bool useDifferentialBeam = true;
    vector3r_t diffBeamDirection = stationDirection;

    // Calculate station beams
    std::cout << ">>> Calculate station beams" << std::endl;
	std::vector<std::complex<float>> aTermBuffer;
    aTermBuffer.resize(subgrid_size*subgrid_size*4*nr_stations);
    calculateStationBeams(stations, itrfDirections, stationDirection, tileDirection, useDifferentialBeam, diffBeamDirection, subgrid_size, aTermBuffer, currentTime, centreFrequency);

    // Store aterm
    std::string aterms_filename("station-beams.fits");
    StoreATermsReal(aterms_filename, aTermBuffer.data(), nr_stations, subgrid_size, subgrid_size);
}
