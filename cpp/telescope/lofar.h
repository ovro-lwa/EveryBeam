// LOFARTelescope.h: Base class for computing the response for the LOFAR
// telescope.
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

#ifndef EVERYBEAM_TELESCOPE_LOFAR_H_
#define EVERYBEAM_TELESCOPE_LOFAR_H_

#include "../station.h"
#include "../elementresponse.h"
#include "telescope.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <memory>

namespace everybeam {

namespace griddedresponse {
class LOFARGrid;
class GriddedResponse;
}  // namespace griddedresponse

namespace telescope {

//! LOFAR telescope class
class LOFAR final : public Telescope {
  friend class griddedresponse::LOFARGrid;

 public:
  /**
   * @brief Construct a new LOFAR object
   *
   * @param ms MeasurementSet
   * @param model Element Response model
   * @param options telescope options
   */
  LOFAR(casacore::MeasurementSet &ms, const ElementResponseModel model,
        const Options &options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) override;

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

 private:
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

  Station::Ptr ReadStation(const casacore::MeasurementSet &ms,
                           const std::size_t id,
                           const ElementResponseModel model) const;
  struct MSProperties {
    double subband_freq;
    casacore::MDirection delay_dir, tile_beam_dir, preapplied_beam_dir;
  };

  std::vector<Station::Ptr> stations_;
  MSProperties ms_properties_;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_LOFAR_H_