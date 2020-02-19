#include <iostream>
#include <complex>
#include <vector>

#include "beam-helper.h"

void calculateElementBeams(
    LOFAR::StationResponse::Station::Ptr& station,
    std::vector<vector2r_t>& thetaPhiDirections,
    size_t nr_antennas,
    unsigned int subgrid_size,
    double frequency,
    std::vector<std::complex<float>>& buffer)
{
    typedef std::complex<float> Data[nr_antennas][subgrid_size][subgrid_size][4];
    Data* data_ptr = (Data *) buffer.data();

    auto elementResponse = station->get_element_response();

    #pragma omp parallel for
    for (size_t a = 0; a < nr_antennas; a++) {
        for (unsigned y = 0; y < subgrid_size; y++) {
            for (unsigned x = 0; x < subgrid_size; x++) {
                // Get theta, phi
                auto direction_thetaphi = thetaPhiDirections[y * subgrid_size + x];
                double theta = direction_thetaphi[0];
                double phi = direction_thetaphi[1];

                // Compute gain
                std::complex<double> gainMatrix[2][2] = {0.0};
                if (std::isfinite(theta) && std::isfinite(phi)) {
                    elementResponse->response(a, frequency, theta, phi, gainMatrix);
                }

                // Store gain
                std::complex<float>* antBufferPtr = (*data_ptr)[a][y][x];
                antBufferPtr[0] = gainMatrix[0][0];
                antBufferPtr[1] = gainMatrix[0][1];
                antBufferPtr[2] = gainMatrix[1][0];
                antBufferPtr[3] = gainMatrix[1][1];
            }
        }
    }
}

void calculateElementBeams(
    LOFAR::StationResponse::Station::Ptr& station,
    std::vector<vector3r_t>& itrfDirections,
    size_t nr_antennas,
    unsigned int subgrid_size,
    double time,
    double frequency,
    std::vector<std::complex<float>>& buffer)
{
    typedef std::complex<float> Data[nr_antennas][subgrid_size][subgrid_size][4];
    Data* data_ptr = (Data *) buffer.data();

    #pragma omp parallel for
    for (size_t a = 0; a < nr_antennas; a++) {
        for (unsigned y = 0; y < subgrid_size; y++) {
            for (unsigned x = 0; x < subgrid_size; x++) {
                // Get direction
                auto direction = itrfDirections[y * subgrid_size + x];

                // Compute gain
                matrix22c_t gainMatrix = { 0.0 };
                if (std::isfinite(direction[0])) {
                    gainMatrix = station->elementResponse(time, frequency, direction, a, true);
                }

                // Store gain
                std::complex<float>* antBufferPtr = (*data_ptr)[a][y][x];
                antBufferPtr[0] = gainMatrix[0][0];
                antBufferPtr[1] = gainMatrix[0][1];
                antBufferPtr[2] = gainMatrix[1][0];
                antBufferPtr[3] = gainMatrix[1][1];
            }
        }
    }
}

void run(
    LOFAR::StationResponse::ElementResponseModel elementResponseModel,
    double frequency,
    std::string& input_filename,
    std::string& output_filename)
{
    // Open measurement set
    std::cout << ">> Opening measurementset: " << input_filename << std::endl;
    casacore::MeasurementSet ms(input_filename);

    // Print frequency
    std::clog << "Frequency: " << frequency * 1e-6 << " Mhz" << std::endl;

    // Set number of stations to 1
    size_t nr_stations = 1;

    // Read number of timesteps
    size_t nr_timesteps = ms.nrow();
    std::clog << "Number of timesteps: " << nr_timesteps << std::endl;


    // Read observation time
    casacore::ScalarColumn<double> timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
    double currentTime = timeColumn(nr_timesteps / 2);

    // Read station
    size_t field_id = 0;
    size_t station_id = 0;
    auto station = readStation(ms, station_id, elementResponseModel);
    auto field_name = GetFieldName(ms, field_id);
    auto station_name = GetStationName(ms, station_id);
    auto nr_antennas = GetNrAntennas(ms, field_id);
    std::cout << "field: " << field_name << std::endl;
    std::cout << "station: " << station_name << std::endl;
    std::cout << "nr_antennas: " << nr_antennas << std::endl;

    // Compute RA and DEC of zenith at currentTime
    double zenithRA, zenithDec;
    GetRaDecZenith(station->position(), currentTime, zenithRA, zenithDec);
    std::clog << "RA:  " << zenithRA << std::endl;
    std::clog << "DEC: " << zenithDec << std::endl;

    // Imaging parameters
    float image_size = M_PI; // in radians
    size_t subgrid_size = 32; // in pixels

    // Compute theta, phi directions
    std::cout << ">>> Computing theta, phi directions to evaluate beam" << std::endl;
    std::vector<vector2r_t> thetaPhiDirections(subgrid_size * subgrid_size);
    GetThetaPhiDirectionsZenith(thetaPhiDirections.data(), subgrid_size);

    // Compute itrfs directions
    std::cout << ">>> Computing theta, phi directions to evaluate beam" << std::endl;
    std::vector<vector3r_t> itrfDirections(subgrid_size * subgrid_size);
    GetITRFDirections(itrfDirections.data(), subgrid_size, image_size, currentTime, zenithRA, zenithDec);

    // Compute element beams from theta, phi
    std::cout << ">>> Computing element beams" << std::endl;
    std::vector<std::complex<float>> beam_thetaphi(subgrid_size*subgrid_size*4*nr_antennas);
    calculateElementBeams(station, thetaPhiDirections, nr_antennas, subgrid_size, frequency, beam_thetaphi);

    // Compute element beams from itrf coordinates
    #if 0
    // TODO: the Station::elementResponse method does not work properly
    std::vector<std::complex<float>> beam_itrf(subgrid_size*subgrid_size*4*nr_antennas);
    calculateElementBeams(station, itrfDirections, nr_antennas, subgrid_size, currentTime, frequency, beam_itrf);
    #endif

    // Store aterm
    std::string aterms_filename(output_filename);
    StoreBeam(aterms_filename, beam_thetaphi.data(), nr_antennas, subgrid_size, subgrid_size);
}
