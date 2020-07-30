// mwa.h: Class for MWA telescopes.
// Inherits from Telescope class.
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

#ifndef EVERYBEAM_TELESCOPE_MWA_H_
#define EVERYBEAM_TELESCOPE_MWA_H_

#include "telescope.h"

#include <casacore/measures/Measures/MPosition.h>

namespace everybeam {

namespace griddedresponse {
class MWAGrid;
class GriddedResponse;
}  // namespace griddedresponse

namespace telescope {

class MWA final : public Telescope {
  friend class griddedresponse::MWAGrid;

 public:
  /**
   * @brief Construct a new MWA object
   *
   * @param ms MeasurementSet
   * @param model Element Response model
   * @param options telescope options
   */
  MWA(casacore::MeasurementSet &ms, const Options &options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) override;

 private:
  struct MSProperties {
    double delays[16];
    casacore::MPosition array_position;
  };
  MSProperties ms_properties_;
};

}  // namespace telescope
}  // namespace everybeam
#endif  // EVERYBEAM_TELESCOPE_MWA_H_