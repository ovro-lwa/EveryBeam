// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DEMO_ELEMENT_BEAM_COMMON_H_
#define DEMO_ELEMENT_BEAM_COMMON_H_

#include <iostream>
#include <complex>
#include <vector>

#include <xtensor/xadapt.hpp>

#include <aocommon/matrix2x2.h>
#include <aocommon/parallelfor.h>
#include <aocommon/system.h>

#include "../cpp/elementresponse.h"
#include "../cpp/options.h"
#include "../cpp/msreadutils.h"
#include "../cpp/station.h"

#include "beam-helper.h"

void calculateElementBeams(
    const everybeam::ElementResponse& elementResponse,
    std::vector<everybeam::vector2r_t>& thetaPhiDirections, size_t nr_antennas,
    unsigned int subgrid_size, double frequency,
    std::vector<std::complex<float>>& buffer) {
  const std::array<size_t, 4> shape{nr_antennas, subgrid_size, subgrid_size, 4};
  auto data = xt::adapt(buffer, shape);

  aocommon::ParallelFor<size_t> loop(aocommon::system::ProcessorCount());
  loop.Run(0, nr_antennas, [&, nr_antennas, subgrid_size](size_t a) {
    for (unsigned y = 0; y < subgrid_size; y++) {
      for (unsigned x = 0; x < subgrid_size; x++) {
        // Get theta, phi
        const everybeam::vector2r_t& direction_thetaphi =
            thetaPhiDirections[y * subgrid_size + x];
        double theta = direction_thetaphi[0];
        double phi = direction_thetaphi[1];

        // Compute gain
        // std::complex<double> gainMatrix[2][2] = {0.0};
        aocommon::MC2x2 gainMatrix;
        if (std::isfinite(theta) && std::isfinite(phi)) {
          gainMatrix = elementResponse.Response(a, frequency, theta, phi);
        }

        // Store gain
        std::complex<float>* antBufferPtr = &data(a, y, x, 0);
        gainMatrix.AssignTo(antBufferPtr);
      }
    }
  });
}

void calculateElementBeams(const everybeam::Station& station,
                           std::vector<everybeam::vector3r_t>& itrfDirections,
                           size_t nr_antennas, unsigned int subgrid_size,
                           double time, double frequency,
                           std::vector<std::complex<float>>& buffer) {
  const std::array<size_t, 4> shape{nr_antennas, subgrid_size, subgrid_size, 4};
  auto data = xt::adapt(buffer, shape);

  aocommon::ParallelFor<size_t> loop(aocommon::system::ProcessorCount());
  loop.Run(0, nr_antennas, [&](size_t a) {
    for (unsigned y = 0; y < subgrid_size; y++) {
      for (unsigned x = 0; x < subgrid_size; x++) {
        // Get direction
        auto direction = itrfDirections[y * subgrid_size + x];

        // Compute gain
        aocommon::MC2x2 gainMatrix(0., 0., 0., 0.);
        if (std::isfinite(direction[0])) {
          gainMatrix = station.ComputeElementResponse(
              time, frequency, direction, a, true, false);
        }

        // Store gain
        std::complex<float>* antBufferPtr = &data(a, y, x, 0);
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
  everybeam::Options options;
  options.element_response_model = elementResponseModel;
  auto station = everybeam::ReadSingleStation(ms, station_id, options);
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
  std::vector<everybeam::vector2r_t> thetaPhiDirections(subgrid_size *
                                                        subgrid_size);
  GetThetaPhiDirectionsZenith(thetaPhiDirections.data(), subgrid_size);
  std::vector<std::complex<float>> beam_thetaphi(subgrid_size * subgrid_size *
                                                 4 * nr_antennas);
  calculateElementBeams(*station->GetElementResponse(), thetaPhiDirections,
                        nr_antennas, subgrid_size, frequency, beam_thetaphi);

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

#endif
