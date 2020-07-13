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

#include "telescope.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <memory>

namespace everybeam {

namespace gridded_response {
class LOFARGrid;
class GriddedResponse;
}  // namespace gridded_response

namespace telescope {

//! LOFAR telescope class
class LOFAR final : public Telescope {
  friend class gridded_response::LOFARGrid;

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

  std::unique_ptr<gridded_response::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) override;

 private:
  Station::Ptr ReadStation(const casacore::MeasurementSet &ms,
                           const std::size_t id,
                           const ElementResponseModel model) const override;
  struct MSProperties {
    double subband_freq;
    casacore::MDirection delay_dir, tile_beam_dir;
  };
  MSProperties ms_properties_;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_LOFAR_H_