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
            const aocommon::CoordinateSystem& coordinate_system)
      : PhasedArrayGrid(telescope_ptr, coordinate_system){};
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_LOFARGRID_H_
