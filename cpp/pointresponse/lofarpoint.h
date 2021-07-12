// lofarpoint.h: Class for computing the LOFAR beam response at given point
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_POINTRESPONSE_LOFARPOINT_H_
#define EVERYBEAM_POINTRESPONSE_LOFARPOINT_H_

#include "phasedarraypoint.h"
#include "../telescope/lofar.h"

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
      : PhasedArrayPoint(telescope_ptr, time) {
    // Extract LOFAR specific options from ms_properties_ and telescope::Options
    // NOTE: same code as in LOFARGrid constructor
    const telescope::LOFAR& lofartelescope =
        dynamic_cast<const telescope::LOFAR&>(*telescope_);
    delay_dir_ = lofartelescope.ms_properties_.delay_dir;
    tile_beam_dir_ = lofartelescope.ms_properties_.tile_beam_dir;
    preapplied_beam_dir_ = lofartelescope.ms_properties_.preapplied_beam_dir;
    preapplied_correction_mode_ =
        lofartelescope.ms_properties_.preapplied_correction_mode;
    subband_frequency_ = lofartelescope.ms_properties_.subband_freq;
    use_channel_frequency_ = lofartelescope.GetOptions().use_channel_frequency;
  }
};
}  // namespace pointresponse
}  // namespace everybeam
#endif  // EVERYBEAM_POINTRESPONSE_LOFARPOINT_H_
