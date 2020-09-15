// lofargrid.h: Class for computing the LOFAR (gridded) response.
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

#include "phasedarraygrid.h"
#include <iostream>
#include <aocommon/matrix2x2.h>
#include <complex>
#include <limits>

namespace everybeam {
namespace griddedresponse {

/**
 * @brief Class for computing the LOFAR gridded response
 *
 */
class LOFARGrid final : public PhasedArrayGrid {
 public:
  /**
   * @brief Construct a new LOFARGrid object
   *
   * @param telescope_ptr Pointer to telescope::LOFAR object
   * @param coordinate_system CoordinateSystem struct
   */
  LOFARGrid(telescope::Telescope* telescope_ptr,
            const coords::CoordinateSystem& coordinate_system)
      : PhasedArrayGrid(telescope_ptr, coordinate_system) {
    // Extract LOFAR specific options from ms_properties_ and telescope::Options
    const telescope::LOFAR& lofartelescope =
        dynamic_cast<const telescope::LOFAR&>(*telescope_);
    delay_dir_ = lofartelescope.ms_properties_.delay_dir;
    tile_beam_dir_ = lofartelescope.ms_properties_.tile_beam_dir;
    preapplied_beam_dir_ = lofartelescope.ms_properties_.preapplied_beam_dir;
    subband_frequency_ = lofartelescope.ms_properties_.subband_freq;
    use_channel_frequency_ = lofartelescope.GetOptions().use_channel_frequency;
  };
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_LOFARGRID_H_
