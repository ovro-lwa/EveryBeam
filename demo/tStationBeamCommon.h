// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DEMO_STATION_BEAM_COMMON_H_
#define DEMO_STATION_BEAM_COMMON_H_

#include <iostream>
#include <complex>
#include <vector>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <xtensor/xadapt.hpp>

#include <aocommon/matrix2x2.h>
#include <aocommon/parallelfor.h>
#include <aocommon/system.h>

#include "../cpp/msreadutils.h"
#include "../cpp/options.h"
#include "../cpp/station.h"

#include "beam-helper.h"

void calculateStationBeams(
    std::vector<std::shared_ptr<everybeam::Station>>& stations,
    std::vector<everybeam::vector3r_t>& itrfDirections,
    everybeam::vector3r_t stationDirection, everybeam::vector3r_t tileDirection,
    unsigned int subgrid_size, std::vector<std::complex<float>>& buffer,
    double time, double frequency) {
  const size_t n_stations = stations.size();
  const std::array<size_t, 4> shape{n_stations, subgrid_size, subgrid_size, 4};
  auto data = xt::adapt(buffer, shape);
  aocommon::ParallelFor<size_t> loop(aocommon::system::ProcessorCount());
  loop.Run(0, stations.size(), [&, subgrid_size, n_stations](size_t s) {
    for (unsigned y = 0; y < subgrid_size; y++) {
      for (unsigned x = 0; x < subgrid_size; x++) {
        auto direction = itrfDirections[y * subgrid_size + x];
        auto freq_beamformer = frequency;
        aocommon::MC2x2 gainMatrix =
            stations[s]->Response(time, frequency, direction, freq_beamformer,
                                  stationDirection, tileDirection);

        std::complex<float>* antBufferPtr = &data(s, y, x, 0);
        gainMatrix.AssignTo(antBufferPtr);
      }
    }
  });
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
  std::vector<std::shared_ptr<everybeam::Station>> stations;
  stations.resize(nr_stations);
  everybeam::Options options;
  options.element_response_model = elementResponseModel;
  everybeam::ReadAllStations(ms, stations.begin(), options);

  // Imaging parameters
  float image_size = 0.5;    // in radians
  size_t subgrid_size = 32;  // in pixels

  // Evaluate beam directions in ITRF coordinates
  std::cout << ">>> Computing directions to evaluate beam" << std::endl;
  std::vector<everybeam::vector3r_t> itrfDirections(subgrid_size *
                                                    subgrid_size);
  GetITRFDirections(itrfDirections.data(), subgrid_size, image_size,
                    currentTime, phaseCentreRA, phaseCentreDec);

  // Set station beam direction to centre of field
  everybeam::vector3r_t stationDirection =
      itrfDirections[(subgrid_size / 2) * subgrid_size + (subgrid_size / 2)];

  // Set tile beam direction equal to station direction
  everybeam::vector3r_t tileDirection = stationDirection;

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

#endif
