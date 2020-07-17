// Telescope.h: Base class for computing the Telescope response.
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

#ifndef EVERYBEAM_TELESCOPE_TELESCOPE_H_
#define EVERYBEAM_TELESCOPE_TELESCOPE_H_

#include "./../station.h"
#include "./../options.h"
#include "./../element_response.h"

#include <vector>
#include <memory>
#include <cassert>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>

namespace everybeam {

namespace gridded_response {
class GriddedResponse;
}

namespace coords {
struct CoordinateSystem;
}

namespace telescope {

/**
 * @brief Telescope class, forming the base class for specific telescopes.
 *
 */
class Telescope {
 public:
  /**
   * @brief Return the gridded response object
   *
   * @param coordinate_system Coordinate system struct
   * @return GriddedResponse::Ptr
   */
  virtual std::unique_ptr<gridded_response::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) = 0;

  /**
   * @brief Get station by index
   *
   * @param station_id Station index to retrieve
   * @return Station::Ptr
   */
  Station::Ptr GetStation(std::size_t station_idx) const {
    // Assert only in DEBUG mode
    assert(station_idx < nstations_);
    return stations_[station_idx];
  }

  std::size_t GetNrStations() const { return nstations_; };

  Options GetOptions() const { return options_; };

 protected:
  /**
   * @brief Construct a new Telescope object
   *
   * @param ms MeasurementSet
   * @param model ElementResponse model
   * @param options telescope options
   */
  Telescope(casacore::MeasurementSet &ms, const ElementResponseModel model,
            const Options &options)
      : nstations_(ms.antenna().nrow()), options_(options) {
    stations_.resize(nstations_);
  };

  /**
   * @brief Read stations into vector
   *
   * @param out_it std::vector of stations
   * @param ms measurement set
   * @param model model
   */
  void ReadAllStations(const casacore::MeasurementSet &ms,
                       const ElementResponseModel model) {
    casacore::ROMSAntennaColumns antenna(ms.antenna());

    for (std::size_t i = 0; i < antenna.nrow(); ++i) {
      stations_[i] = ReadStation(ms, i, model);
    }
  };

  std::size_t nstations_;
  Options options_;
  std::vector<Station::Ptr> stations_;

 private:
  // Virtual method for reading a single telescope station
  virtual Station::Ptr ReadStation(const casacore::MeasurementSet &ms,
                                   const std::size_t id,
                                   const ElementResponseModel model) const = 0;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_TELESCOPE_H_