// dishgrid.h: Class for computing the circular symmetric (gridded) response.
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

#ifndef EVERYBEAM_GRIDDEDRESPONSE_DISHGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_DISHGRID_H_

#include "griddedresponse.h"

namespace everybeam {
namespace griddedresponse {

class DishGrid final : public GriddedResponse {
 public:
  DishGrid(telescope::Telescope* telescope_ptr,
           const coords::CoordinateSystem coordinate_system)
      : GriddedResponse(telescope_ptr, coordinate_system){};

  void CalculateStation(std::complex<float>* buffer, double time,
                        double frequency, size_t station_idx,
                        size_t field_id) override;

  void CalculateAllStations(std::complex<float>* buffer, double time,
                            double frequency, size_t field_id) override;
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_DISHGRID_H_