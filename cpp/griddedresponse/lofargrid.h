// lofargrid.h: Class for computing the LOFAR (gridded) response.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_LOFARGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_LOFARGRID_H_

#include "phasedarraygrid.h"
#include "../telescope/lofar.h"

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
  LOFARGrid(const telescope::Telescope* telescope_ptr,
            const coords::CoordinateSystem& coordinate_system)
      : PhasedArrayGrid(telescope_ptr, coordinate_system) {
    // Extract LOFAR specific options from ms_properties_ and telescope::Options
    const telescope::LOFAR& lofartelescope =
        dynamic_cast<const telescope::LOFAR&>(*telescope_);
    delay_dir_ = lofartelescope.ms_properties_.delay_dir;
    tile_beam_dir_ = lofartelescope.ms_properties_.tile_beam_dir;
    preapplied_beam_dir_ = lofartelescope.ms_properties_.preapplied_beam_dir;
    preapplied_correction_mode_ =
        lofartelescope.ms_properties_.preapplied_correction_mode;
    subband_frequency_ = lofartelescope.ms_properties_.subband_freq;
    use_channel_frequency_ = lofartelescope.GetOptions().use_channel_frequency;
  };
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_LOFARGRID_H_
