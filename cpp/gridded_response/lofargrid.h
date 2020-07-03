// LOFARGrid.h: Class for computing the LOFAR (gridded) response.
//
// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the EveryBeam software suite.
// The EveryBeam software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The EveryBeam software suite is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the EveryBeam software suite. If not, see
// <http://www.gnu.org/licenses/>.
//
// $Id$

#ifndef EVERYBEAM_GRIDDEDRESPONSE_LOFARGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_LOFARGRID_H_

#include "griddedresponse.h"
#include <iostream>
#include <complex>

namespace everybeam {
namespace gridded_response {

/**
 * @brief Class for computing the LOFAR gridded response
 *
 */
class LOFARGrid final : public GriddedResponse {
 public:
  typedef std::unique_ptr<LOFARGrid> Ptr;

  /**
   * @brief Construct a new LOFARGrid object
   *
   * @param telescope_ptr Pointer to telescope::LOFAR object
   * @param coordinateSystem CoordinateSystem struct
   */
  LOFARGrid(telescope::Telescope* telescope_ptr,
            const CoordinateSystem& coordinateSystem)
      : GriddedResponse(telescope_ptr, coordinateSystem){};

  // TODO: remove this placeholders
  void CalculateStation() { std::cout << "One is good" << std::endl; };

  // TODO: remove this placeholder
  void CalculateStation(std::size_t station_id) {
    auto station_tmp = _telescope->GetStation(station_id);
    std::cout << "Station name for index " << station_id << " is "
              << station_tmp->name() << std::endl;
  };

  // TODO: remove this placeholder
  void CalculateAllStations() {
    std::size_t val = 0;
    for (std::size_t station_id = 0; station_id < _telescope->stations.size();
         ++station_id) {
      val++;
      // Repeated call to CalculateStation?
      auto station_tmp = _telescope->GetStation(station_id);
    };
    std::cout << "I just read " << val << " stations" << std::endl;
  };

  // Geared towards implementation
  /**
   * @brief Compute the Beam response for a single station
   *
   * @param buffer Output buffer
   * @param stationIdx Station index, must be smaller than number of stations
   * in the Telescope
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency (Hz)
   */
  void CalculateStation(std::complex<float>* buffer, size_t stationIdx,
                        double time, double freq) override{};
  // Repeated call of calculate single?
  /**
   * @brief Compute the Beam response for all stations in a Telescope
   *
   * @param buffer Output buffer
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency (Hz)
   */
  void CalculateAllStations(std::complex<float>* buffer, double time,
                            double freq) override{};
};
}  // namespace gridded_response
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_LOFARGRID_H_