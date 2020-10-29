// oskar.h: Base class for computing the response for the OSKAR
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

#ifndef EVERYBEAM_TELESCOPE_OSKAR_H_
#define EVERYBEAM_TELESCOPE_OSKAR_H_

#include "../station.h"
#include "../elementresponse.h"
#include "phasedarray.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <memory>

namespace everybeam {

namespace griddedresponse {
class OSKARGrid;
class GriddedResponse;
}  // namespace griddedresponse

namespace telescope {

//! OSKAR telescope class
class OSKAR final : public PhasedArray {
  friend class griddedresponse::OSKARGrid;

 public:
  /**
   * @brief Construct a new OSKAR object
   *
   * @param ms MeasurementSet
   * @param model Element Response model
   * @param options telescope options
   */
  OSKAR(casacore::MeasurementSet &ms, const Options &options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) override;

  //! Get the tile beam direction, equal to delay direction for OSKAR!
  virtual casacore::MDirection GetTileBeamDirection() const final override {
    std::cout << "OSKAR has no tile. tile_beam_dir is equal to the delay_dir."
              << std::endl;
    return ms_properties_.tile_beam_dir;
  };

  //! Get the preapplied beam direction, equal to delay direction for OSKAR!
  virtual casacore::MDirection GetPreappliedBeamDirection()
      const final override {
    std::cout << "OSKAR has no preapplied beam direction (yet). "
                 "preapplied_beam_dir is equal to the delay_dir."
              << std::endl;
    return ms_properties_.preapplied_beam_dir;
  };

 private:
  Station::Ptr ReadStation(const casacore::MeasurementSet &ms,
                           const std::size_t id,
                           const ElementResponseModel model) const override;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_OSKAR_H_
