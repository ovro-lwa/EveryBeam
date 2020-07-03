#include <boost/test/unit_test.hpp>

#include "./../load.h"
#include "./../options.h"
#include "./../ElementResponse.h"

#include "config.h"

using namespace everybeam;

BOOST_AUTO_TEST_CASE(load_lofar) {
  ElementResponseModel response_model = ElementResponseModel::Hamaker;
  Options options;
  casacore::MeasurementSet ms(LOFAR_MOCK_MS);

  // Load LOFAR Telescope
  auto telescope = Load(ms, response_model, options);

  // Assert if we indeed have a LOFAR pointer
  BOOST_CHECK(nullptr != dynamic_cast<telescope::LOFAR*>(telescope.get()));

  // Assert if correct number of stations
  std::size_t nstations = 70;
  BOOST_CHECK_EQUAL(telescope->stations.size(), nstations);

  // Assert if GetStation(stationd_id) behaves properly
  BOOST_CHECK_EQUAL(telescope->GetStation(0)->name(), "CS001HBA0");

  // Get gridded response
  CoordinateSystem coord_system;
  auto grrp = telescope->GetGriddedResponse(coord_system);
  BOOST_CHECK(nullptr !=
              dynamic_cast<gridded_response::LOFARGrid*>(grrp.get()));

  // TODO: add test
  //   grrp->CalculateStation(0);
  //   grrp->CalculateStation();
  //   grrp->CalculateAllStations();
}