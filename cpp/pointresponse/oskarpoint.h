// oskarpoint.h: Class for computing the OSKAR response at given point
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_POINTRESPONSE_OSKARPOINT_H_
#define EVERYBEAM_POINTRESPONSE_OSKARPOINT_H_

#include "phasedarraypoint.h"
#include "../telescope/oskar.h"

namespace everybeam {
namespace pointresponse {

/**
 * @brief Class for computing the OSKAR SKALA 4.0 response for given ra/dec
 *
 */
class OSKARPoint final : public PhasedArrayPoint {
 public:
  /**
   * @brief Construct a new OSKARPoint object
   *
   * @param telescope_ptr Pointer to telescope::OSKAR object
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s)
   */
  OSKARPoint(const telescope::Telescope* telescope_ptr, double time)
      : PhasedArrayPoint(telescope_ptr, time) {
    // NOTE: for OSKAR, it always holds that:
    use_channel_frequency_ = true;
    // i.e. never use the subband frequency (which is explicitly set to 0):
    subband_frequency_ = 0.0;
  }
};
}  // namespace pointresponse
}  // namespace everybeam
#endif  // EVERYBEAM_POINTRESPONSE_OSKARPOINT_H_