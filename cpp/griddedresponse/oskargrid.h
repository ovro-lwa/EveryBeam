// oskargrid.h: Class for computing the OSKAR (gridded) response.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_OSKARGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_OSKARGRID_H_

#include "phasedarraygrid.h"
#include "../telescope/oskar.h"

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
  OSKARGrid(const telescope::Telescope* telescope_ptr,
            const aocommon::CoordinateSystem& coordinate_system)
      : PhasedArrayGrid(telescope_ptr, coordinate_system) {
    // NOTE: for OSKAR, it always holds that:
    use_channel_frequency_ = true;
    // i.e. never use the subband frequency (which is explicitly set to 0):
    subband_frequency_ = 0.0;
  };
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_OSKARGRID_H_
