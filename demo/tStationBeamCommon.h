#include <iostream>
#include <complex>
#include <vector>

#include "beam-helper.h"

void calculateStationBeams(std::vector<everybeam::Station::Ptr>& stations,
                           std::vector<vector3r_t>& itrfDirections,
                           vector3r_t stationDirection,
                           vector3r_t tileDirection, unsigned int subgrid_size,
                           std::vector<std::complex<float>>& buffer,
                           double time, double frequency) {
  typedef std::complex<float> Data[stations.size()][subgrid_size][subgrid_size]
                                  [4];
  Data* data_ptr = (Data*)buffer.data();

#pragma omp parallel for
  for (size_t s = 0; s < stations.size(); s++) {
    for (unsigned y = 0; y < subgrid_size; y++) {
      for (unsigned x = 0; x < subgrid_size; x++) {
        auto direction = itrfDirections[y * subgrid_size + x];
        auto freq_beamformer = frequency;
        matrix22c_t gainMatrix =
            stations[s]->Response(time, frequency, direction, freq_beamformer,
                                  stationDirection, tileDirection);

        std::complex<float>* antBufferPtr = (*data_ptr)[s][y][x];

        matrix22c_t stationGain = gainMatrix;
        antBufferPtr[0] = stationGain[0][0];
        antBufferPtr[1] = stationGain[0][1];
        antBufferPtr[2] = stationGain[1][0];
        antBufferPtr[3] = stationGain[1][1];
      }
    }
  }
}

void run(everybeam::ElementResponseModel elementResponseModel, double frequency,
         std::string& input_filename, std::string& output_filename) {
  // Open measurement set
  std::cout << ">> Opening measurementset: " << input_filename << std::endl;
  casacore::MeasurementSet ms(input_filename);

  // Print frequency
  std::clog << "Frequency: " << frequency * 1e-6 << " Mhz" << std::endl;

  // Read number of stations
  size_t nr_stations = ms.antenna().nrow();
  std::clog << "Number of stations: " << nr_stations << std::endl;

  // Read number of timesteps
  size_t nr_timesteps = ms.nrow();
  std::clog << "Number of timesteps: " << nr_timesteps << std::endl;

  // Read field id
  casacore::ROScalarColumn<int> fieldIdColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
  size_t fieldId = fieldIdColumn.getColumn()[0];
  std::clog << "Field ID: " << fieldId << std::endl;

  // Read phase centre info
  double phaseCentreRA, phaseCentreDec;
  GetPhaseCentreInfo(ms, fieldId, phaseCentreRA, phaseCentreDec);
  std::clog << "RA:  " << phaseCentreRA << std::endl;
  std::clog << "DEC: " << phaseCentreDec << std::endl;

  // Read observation time
  casacore::ScalarColumn<double> timeColumn(
      ms, ms.columnName(casacore::MSMainEnums::TIME));
  double currentTime = timeColumn(nr_timesteps / 2);

  // Read stations
  std::vector<everybeam::Station::Ptr> stations;
  stations.resize(nr_stations);
  ReadStations(ms, stations.begin(), elementResponseModel);

  // Imaging parameters
  float image_size = 0.5;    // in radians
  size_t subgrid_size = 32;  // in pixels

  // Evaluate beam directions in ITRF coordinates
  std::cout << ">>> Computing directions to evaluate beam" << std::endl;
  std::vector<vector3r_t> itrfDirections(subgrid_size * subgrid_size);
  GetITRFDirections(itrfDirections.data(), subgrid_size, image_size,
                    currentTime, phaseCentreRA, phaseCentreDec);

  // Set station beam direction to centre of field
  vector3r_t stationDirection =
      itrfDirections[(subgrid_size / 2) * subgrid_size + (subgrid_size / 2)];

  // Set tile beam direction equal to station direction
  vector3r_t tileDirection = stationDirection;

  // Compute station beams
  std::cout << ">>> Computing station beams" << std::endl;
  std::vector<std::complex<float>> beams;
  beams.resize(subgrid_size * subgrid_size * 4 * nr_stations);
  calculateStationBeams(stations, itrfDirections, stationDirection,
                        tileDirection, subgrid_size, beams, currentTime,
                        frequency);

  // Store aterm
  std::cout << ">>> Writing beam images to: " << output_filename << std::endl;
  StoreBeam(output_filename, beams.data(), nr_stations, subgrid_size,
            subgrid_size);
}
