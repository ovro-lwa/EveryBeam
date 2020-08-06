#include <iostream>
#include <complex>
#include <vector>

#include "beam-helper.h"

void calculateElementBeams(everybeam::Station::Ptr& station,
                           std::vector<vector2r_t>& thetaPhiDirections,
                           size_t nr_antennas, unsigned int subgrid_size,
                           double frequency,
                           std::vector<std::complex<float>>& buffer) {
  typedef std::complex<float> Data[nr_antennas][subgrid_size][subgrid_size][4];
  Data* data_ptr = (Data*)buffer.data();

  auto elementResponse = station->GetElementResponse();

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
          elementResponse->Response(a, frequency, theta, phi, gainMatrix);
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

void calculateElementBeams(everybeam::Station::Ptr& station,
                           std::vector<vector3r_t>& itrfDirections,
                           size_t nr_antennas, unsigned int subgrid_size,
                           double time, double frequency,
                           std::vector<std::complex<float>>& buffer) {
  typedef std::complex<float> Data[nr_antennas][subgrid_size][subgrid_size][4];
  Data* data_ptr = (Data*)buffer.data();

#pragma omp parallel for
  for (size_t a = 0; a < nr_antennas; a++) {
    for (unsigned y = 0; y < subgrid_size; y++) {
      for (unsigned x = 0; x < subgrid_size; x++) {
        // Get direction
        auto direction = itrfDirections[y * subgrid_size + x];

        // Compute gain
        matrix22c_t gainMatrix = {0.0};
        if (std::isfinite(direction[0])) {
          gainMatrix = station->ComputeElementResponse(time, frequency,
                                                       direction, a, true);
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

void run(everybeam::ElementResponseModel elementResponseModel, double frequency,
         std::string& input_filename, std::string& output_filename) {
  // Open measurement set
  std::cout << ">> Opening measurementset: " << input_filename << std::endl;
  casacore::MeasurementSet ms(input_filename);

  // Print frequency
  std::clog << "Frequency: " << frequency * 1e-6 << " Mhz" << std::endl;

  // Read number of timesteps
  size_t nr_timesteps = ms.nrow();
  std::clog << "Number of timesteps: " << nr_timesteps << std::endl;

  // Read observation time
  casacore::ScalarColumn<double> timeColumn(
      ms, ms.columnName(casacore::MSMainEnums::TIME));
  double currentTime = timeColumn(nr_timesteps / 2);

  // Read station
  size_t field_id = 0;
  size_t station_id = 0;
  auto station = ReadLofarStation(ms, station_id, elementResponseModel);
  auto field_name = GetFieldName(ms, field_id);
  auto station_name = GetStationName(ms, station_id);
  auto nr_antennas = GetNrAntennas(ms, field_id);
  std::cout << "field: " << field_name << std::endl;
  std::cout << "station: " << station_name << std::endl;
  std::cout << "nr_antennas: " << nr_antennas << std::endl;

  // Compute RA and DEC of zenith at currentTime
  double zenithRA, zenithDec;
  GetRaDecZenith(station->GetPosition(), currentTime, zenithRA, zenithDec);
  std::clog << "RA:  " << zenithRA << std::endl;
  std::clog << "DEC: " << zenithDec << std::endl;

  // Imaging parameters
  // float image_size = M_PI;   // in radians
  size_t subgrid_size = 32;  // in pixels

  // Compute element beams from theta, phi
  std::cout << ">>> Computing element beams theta, phi" << std::endl;
  std::vector<vector2r_t> thetaPhiDirections(subgrid_size * subgrid_size);
  GetThetaPhiDirectionsZenith(thetaPhiDirections.data(), subgrid_size);
  std::vector<std::complex<float>> beam_thetaphi(subgrid_size * subgrid_size *
                                                 4 * nr_antennas);
  calculateElementBeams(station, thetaPhiDirections, nr_antennas, subgrid_size,
                        frequency, beam_thetaphi);

// Compute element beams from itrf coordinates
// TODO: the Station::ComputeElementResponse method does not work properly
#if 0
    std::cout << ">>> Computing element beams itrfs" << std::endl;
    std::vector<vector3r_t> itrfDirections(subgrid_size * subgrid_size);
    GetITRFDirections(itrfDirections.data(), subgrid_size, image_size, currentTime, zenithRA, zenithDec);
    std::vector<std::complex<float>> beam_itrf(subgrid_size*subgrid_size*4*nr_antennas);
    calculateElementBeams(station, itrfDirections, nr_antennas, subgrid_size, currentTime, frequency, beam_itrf);
#endif

  // Store aterm
  StoreBeam(output_filename, beam_thetaphi.data(), nr_antennas, subgrid_size,
            subgrid_size);
}
