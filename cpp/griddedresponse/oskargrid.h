// oskargrid.h: Class for computing the OSKAR (gridded) response.
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

#ifndef EVERYBEAM_GRIDDEDRESPONSE_OSKARGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_OSKARGRID_H_

#include "phasedarraygrid.h"

namespace everybeam {
namespace griddedresponse {

/**
 * @brief Class for computing the OSKAR gridded response
 *
 */
class OSKARGrid final : public PhasedArrayGrid {
 public:
  /**
   * @brief Construct a new OSKARGrid object
   *
   * @param telescope_ptr Pointer to telescope::OSKAR object
   * @param coordinate_system CoordinateSystem struct
   */
  OSKARGrid(telescope::Telescope* telescope_ptr,
            const coords::CoordinateSystem& coordinate_system)
      : PhasedArrayGrid(telescope_ptr, coordinate_system) {
    const telescope::OSKAR& oskartelescope =
        dynamic_cast<const telescope::OSKAR&>(*telescope_);

    delay_dir_ = oskartelescope.ms_properties_.delay_dir;
    tile_beam_dir_ = oskartelescope.ms_properties_.delay_dir;
    // NOTE: for OSKAR, it always holds that:
    // use_channel_frequency_ = true (i.e. never use the subband_frequency_)
  };
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_OSKARGRID_H_
