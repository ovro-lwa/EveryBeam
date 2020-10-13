// phasedarray.h: Base class for computing the response for phased array
// telescopes (e.g. LOFAR, SKA-LOW)
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

#ifndef EVERYBEAM_TELESCOPE_PHASEDARRAY_H_
#define EVERYBEAM_TELESCOPE_PHASEDARRAY_H_

#include "../station.h"
#include "../elementresponse.h"
#include "telescope.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <memory>

namespace everybeam {
namespace telescope {

//! PhasedArray telescope class, is parent to OSKAR and LOFAR
class PhasedArray : public Telescope {
 public:
  /**
   * @brief Construct a new PhasedArray object
   *
   * @param ms MeasurementSet
   * @param model Element Response model
   * @param options telescope options
   */
  PhasedArray(casacore::MeasurementSet &ms, const Options &options)
      : Telescope(ms, options) {
    stations_.resize(nstations_);
  };

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

 protected:
  std::vector<Station::Ptr> stations_;

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

  virtual Station::Ptr ReadStation(const casacore::MeasurementSet &ms,
                                   const std::size_t id,
                                   const ElementResponseModel model) const = 0;
};
}  // namespace telescope
}  // namespace everybeam
#endif  // EVERYBEAM_TELESCOPE_PHASEDARRAY_H_