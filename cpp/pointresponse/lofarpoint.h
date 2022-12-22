// lofarpoint.h: Class for computing the LOFAR beam response at given point
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_POINTRESPONSE_LOFARPOINT_H_
#define EVERYBEAM_POINTRESPONSE_LOFARPOINT_H_

#include "phasedarraypoint.h"

namespace everybeam {
namespace pointresponse {

/**
 * @brief Class for computing the LOFAR response for given ra/dec
 *
 */
class LOFARPoint final : public PhasedArrayPoint {
 public:
  /**
   * @param telescope_ptr Pointer to telescope::LOFAR object
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s)
   */
  LOFARPoint(const telescope::Telescope* telescope_ptr, double time)
      : PhasedArrayPoint(telescope_ptr, time) {}
};
}  // namespace pointresponse
}  // namespace everybeam
#endif  // EVERYBEAM_POINTRESPONSE_LOFARPOINT_H_
