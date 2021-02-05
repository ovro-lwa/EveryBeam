// phasedarraypoint.h: class for computing the directional telescope
// responses for OSKAR and LOFAR telescope(s)
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_POINTRESPONSE_PHASEDARRAYPOINT_H_
#define EVERYBEAM_POINTRESPONSE_PHASEDARRAYPOINT_H_

#include "pointresponse.h"
#include "../common/types.h"

#include <aocommon/matrix2x2.h>
#include <casacore/measures/Measures/MDirection.h>

#include <mutex>

namespace everybeam {
namespace pointresponse {

class PhasedArrayPoint : public PointResponse {
 public:
  PhasedArrayPoint(const telescope::Telescope* telescope_ptr, double time);

  void CalculateStation(std::complex<float>* buffer, double ra, double dec,
                        double freq, size_t station_idx,
                        size_t field_id) final override;

 protected:
  casacore::MDirection delay_dir_, tile_beam_dir_, preapplied_beam_dir_;
  vector3r_t station0_, tile0_, dir_itrf_, diff_beam_centre_;

  bool use_channel_frequency_, use_differential_beam_;
  double subband_frequency_;

 private:
  aocommon::MC2x2F inverse_central_gain_;

  /**
   * @brief Method for computing the ITRF-vectors, given a ra, dec position in
   * radians. NOTE: method not thread-safe due to casacore dependencies.
   */
  void SetITRFVectors(double ra, double dec);

  std::mutex mtx_;
};

}  // namespace pointresponse
}  // namespace everybeam

#endif  // EVERYBEAM_POINTRESPONSE_PHASEDARRAYPOINT_H_