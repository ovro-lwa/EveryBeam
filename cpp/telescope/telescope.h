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
#include "./../gridded_response/griddedresponse.h"

#include <vector>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>

namespace everybeam {

struct CoordinateSystem;

namespace telescope {

/**
 * @brief Telescope class, forming the base class for specific telescopes.
 *
 */
class Telescope {
 public:
  typedef std::unique_ptr<Telescope> Ptr;

  // Will be filled by the private read_all_stations() method
  std::vector<Station::Ptr> stations;

  /**
   * @brief Return the gridded response object
   *
   * @param coordinate_system Coordinate system struct
   * @return GriddedResponse::Ptr
   */
  virtual gridded_response::GriddedResponse::Ptr GetGriddedResponse(
      const CoordinateSystem &coordinate_system) = 0;

  /**
   * @brief Get station by index
   *
   * @param station_id Station index to retrieve
   * @return Station::Ptr
   */
  Station::Ptr GetStation(std::size_t station_id) const {
    // TODO: throw exception if idx >_nstations

    return stations[station_id];
  }

 protected:
  /**
   * @brief Construct a new Telescope object
   *
   * @param ms MeasurementSet
   * @param model ElementResponse model
   * @param options telescope options
   */
  Telescope(const casacore::MeasurementSet &ms,
            const ElementResponseModel model, const Options &options)
      : _nstations(ms.antenna().nrow()), _options(options) {
    stations.resize(_nstations);
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
      this->stations[i] = ReadStation(ms, i, model);
    }
  };

  std::size_t _nstations;
  Options _options;

 private:
  // Virtual method for reading a single telescope station
  virtual Station::Ptr ReadStation(const casacore::MeasurementSet &ms,
                                   const std::size_t id,
                                   const ElementResponseModel model) const = 0;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_TELESCOPE_H_