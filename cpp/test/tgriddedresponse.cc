// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "../load.h"
#include "../options.h"
#include "../beammode.h"
#include "../griddedresponse/lofargrid.h"
#include "../pointresponse/lofarpoint.h"
#include "../elementresponse.h"
#include "../coords/coordutils.h"
#include "../station.h"
#include "../common/types.h"
#include "../telescope/lofar.h"
#include "../aterms/atermconfig.h"
#include "../aterms/parsetprovider.h"
#include "../msreadutils.h"

#include "config.h"
#include <complex>
#include <cmath>
#include <iostream>
#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>

using aocommon::CoordinateSystem;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::LOFARGrid;
using everybeam::telescope::LOFAR;
using everybeam::telescope::Telescope;

BOOST_AUTO_TEST_SUITE(gridded_response)

BOOST_AUTO_TEST_CASE(undersampled) {
  everybeam::Options options;
  options.element_response_model = everybeam::ElementResponseModel::kHamaker;

  casacore::MeasurementSet ms(LOFAR_HBA_MOCK_MS);
  std::unique_ptr<Telescope> telescope = Load(ms, options);

  const std::vector<double> times{4929192878.008341};
  constexpr double frequency = 138476562.5;

  constexpr size_t sampling_factor = 8;
  constexpr size_t width = 100;
  constexpr size_t height = 70;
  CoordinateSystem coord_system;
  coord_system.width = width;
  coord_system.height = height;
  coord_system.ra = 2.15374123;
  coord_system.dec = 0.8415521;
  coord_system.dl = 0.05 * M_PI / 180.0;
  coord_system.dm = 0.05 * M_PI / 180.0;
  coord_system.l_shift = 0.0;
  coord_system.m_shift = 0.0;

  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(coord_system);

  std::vector<float> reference_data(grid_response->GetIntegratedBufferSize());
  const std::vector<double> baseline_weights(
      telescope->GetNrStations() * (telescope->GetNrStations() + 1) / 2, 1.0);

  // Use the simple interface to calculate the reference values
  grid_response->IntegratedResponse(everybeam::BeamMode::kFull,
                                    reference_data.data(), times, frequency, 0,
                                    sampling_factor, baseline_weights);

  // Now recalculate in separate steps
  const std::vector<aocommon::HMC4x4> undersampled =
      grid_response->UndersampledIntegratedResponse(
          everybeam::BeamMode::kFull, times, frequency, 0, sampling_factor,
          baseline_weights);

  for (size_t element_index = 0; element_index != 16; ++element_index) {
    std::vector<float> result(width * height);
    GriddedResponse::UpsampleResponse(result.data(), element_index, width,
                                      height, undersampled, sampling_factor);
    for (size_t i = 0; i != result.size(); ++i) {
      BOOST_CHECK_EQUAL(result[i],
                        reference_data[i + element_index * width * height]);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
