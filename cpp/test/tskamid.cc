// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include "config.h"
#include "../load.h"
#include "../options.h"
#include "../griddedresponse/skamidgrid.h"
#include "../pointresponse/skamidpoint.h"

#include <aocommon/imagecoordinates.h>

#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>

using aocommon::CoordinateSystem;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::griddedresponse::SkaMidGrid;
using everybeam::pointresponse::PointResponse;
using everybeam::pointresponse::SkaMidPoint;

namespace everybeam {

namespace {
// Time is currently not relevant
const double kTime = 0.0;
const double kFrequency = 1.3e9;

const double kRa = 0.;
const double kDec = -0.785398;
const std::size_t kWidth = 4;
const std::size_t kHeight = 4;
const double kDl = 0.125 * M_PI / 180;
const double kDm = kDl;
const double kShiftL = 0.;
const double kShiftM = 0.;

CoordinateSystem kCoordSystem = {kWidth, kHeight, kRa,     kDec,
                                 kDl,    kDm,     kShiftL, kShiftM};

}  // namespace

BOOST_AUTO_TEST_SUITE(ska_mid)

BOOST_AUTO_TEST_CASE(compare_responses, *boost::unit_test::tolerance(1e-5)) {
  Options options;
  options.element_response_model = ElementResponseModel::kSkaMidAnalytical;

  casacore::MeasurementSet ms(SKA_MID_MOCK_MS);

  std::unique_ptr<telescope::Telescope> telescope = Load(ms, options);
  std::unique_ptr<GriddedResponse> grid_response =
      telescope->GetGriddedResponse(kCoordSystem);

  // Define buffer and get gridded response (for a single station)
  const size_t n_stations = 1;
  const size_t station_index = 0;

  std::vector<std::complex<float>> antenna_buffer_single(
      grid_response->GetStationBufferSize(n_stations));

  grid_response->Response(everybeam::BeamMode::kFull,
                          antenna_buffer_single.data(), kTime, kFrequency,
                          station_index, 0);

  std::array<std::complex<float>, 4> point_response_buffer;
  std::unique_ptr<PointResponse> point_response =
      telescope->GetPointResponse(kTime);

  // Loop over grid
  for (size_t y = 0; y != kCoordSystem.height; ++y) {
    for (size_t x = 0; x != kCoordSystem.width; ++x) {
      double l;
      double m;
      double ra;
      double dec;
      aocommon::ImageCoordinates::XYToLM(x, y, kCoordSystem.dl, kCoordSystem.dm,
                                         kCoordSystem.width,
                                         kCoordSystem.height, l, m);
      aocommon::ImageCoordinates::LMToRaDec(l, m, kCoordSystem.ra,
                                            kCoordSystem.dec, ra, dec);
      point_response->Response(everybeam::BeamMode::kFull,
                               point_response_buffer.data(), ra, dec,
                               kFrequency, 0, 0);
      const size_t buffer_offset = (y * kCoordSystem.width + x) * 4;
      // Compare gridded response with point response
      BOOST_CHECK_EQUAL_COLLECTIONS(
          point_response_buffer.begin(), point_response_buffer.end(),
          antenna_buffer_single.begin() + buffer_offset,
          antenna_buffer_single.begin() + buffer_offset + 4);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace everybeam
