#include <iostream>
#include <complex>
#include <vector>

#include "beam-helper.h"

void calculateElementBeams(
	LOFAR::StationResponse::Station::Ptr& station,
    size_t nr_antennas,
    std::vector<vector3r_t>& itrfsDirections,
    vector3r_t stationDirection,
    vector3r_t tileDirection,
    bool useDifferentialBeam,
    vector3r_t diffBeamCentre,
    unsigned int subgrid_size,
    std::vector<std::complex<float>>& buffer,
    double time, double frequency)
{
    typedef std::complex<float> Data[nr_antennas][subgrid_size][subgrid_size][4];
    Data* data_ptr = (Data *) buffer.data();

    #pragma omp parallel for
	for (size_t a = 0; a < nr_antennas; a++) {
        for (unsigned y = 0; y < subgrid_size; y++) {
            for (unsigned x = 0; x < subgrid_size; x++) {
                auto direction = itrfsDirections[y * subgrid_size + x];
                auto freq_beamformer = frequency;
                bool rotate = true;
                matrix22c_t gainMatrix =
                    station->elementResponse(
                    time, frequency, direction,
                    a, rotate);

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

    // Read station
    size_t station_id = 0;
    LOFAR::StationResponse::Station::Ptr station = readStation(ms, station_id, elementResponseModel);
    casacore::Table fieldTable = ms.keywordSet().asTable("LOFAR_ANTENNA_FIELD");
    casacore::ROScalarColumn<casacore::String> c_name(fieldTable, "NAME");
    casacore::ROArrayQuantColumn<casacore::Double> c_offset(fieldTable, "ELEMENT_OFFSET", "m");
    const string &name = c_name(station_id);
    casacore::Matrix<casacore::Quantity> aips_offset = c_offset(station_id);
    size_t nr_antennas = aips_offset.ncolumn();
    std::cout << "station " << station_id << ": " << name << std::endl;
    std::cout << "nr_antennas: " << nr_antennas << std::endl;

    // Imaging parameters
    float image_size = 2; // in radians
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

    // Compute station beams
    std::cout << ">>> Computing element beams" << std::endl;
	std::vector<std::complex<float>> aTermBuffer;
    aTermBuffer.resize(subgrid_size*subgrid_size*4*nr_antennas);
    calculateElementBeams(station, nr_antennas, itrfDirections, stationDirection, tileDirection, useDifferentialBeam, diffBeamDirection, subgrid_size, aTermBuffer, currentTime, frequency);

    // Store aterm
    std::string aterms_filename(output_filename);
    StoreBeam(aterms_filename, aTermBuffer.data(), nr_antennas, subgrid_size, subgrid_size);
}
