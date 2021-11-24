// aartfaacpoint.h: Class for computing the point response for
// LOFAR observations in AARTFAAC modus
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_AARTFAACPOINT_H_
#define EVERYBEAM_GRIDDEDRESPONSE_AARTFAACPOINT_H_

#include "../telescope/telescope.h"
#include "phasedarraypoint.h"

namespace everybeam {
namespace pointresponse {
/**
 * @brief Class for computing the AARTFAAC response
 * for given (ra, dec) coordinate. Class exploits the
 * fact that the station response is identical for all stations
 * within an AARTFAAC observation.
 *
 */
class AartfaacPoint final : public PhasedArrayPoint {
 public:
  /**
   * @param telescope_ptr Pointer to telescope::LOFAR object
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s)
   */
  AartfaacPoint(const telescope::Telescope* telescope_ptr, double time)
      : PhasedArrayPoint(telescope_ptr, time) {}

  void ResponseAllStations(BeamMode beam_mode, std::complex<float>* buffer,
                           double ra, double dec, double freq,
                           size_t field_id) override {
    Response(beam_mode, buffer, ra, dec, freq, 0u, field_id);
    // Just repeat nstations times
    for (size_t i = 1; i != telescope_->GetNrStations(); ++i) {
      std::copy_n(buffer, 4, buffer + i * 4);
    }
  }
};
}  // namespace pointresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_AARTFAACPOINT_H_