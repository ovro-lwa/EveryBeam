// phasedarraypoint.h: class for computing the directional telescope
// responses for OSKAR and LOFAR telescope(s)
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_POINTRESPONSE_PHASEDARRAYPOINT_H_
#define EVERYBEAM_POINTRESPONSE_PHASEDARRAYPOINT_H_

#include "pointresponse.h"
#include "../common/types.h"
#include "../correctionmode.h"

#include <aocommon/matrix2x2.h>
#include <casacore/measures/Measures/MDirection.h>

namespace everybeam {
namespace pointresponse {

class PhasedArrayPoint : public PointResponse {
 public:
  PhasedArrayPoint(const telescope::Telescope* telescope_ptr, double time);

  /**
   * @brief Get beam response for a given station at a prescribed ra, dec
   * position.
   * NOTE: the \param ra, \param dec input values are only used if values are
   * different from the cached values. Direction values in cache along with the
   * ITRF directions can be precomputed with UpdateITRFVectors for efficiency.
   * NOTE: CalculateStation complies with the standard
   * threading rules, but does not guarantee thread-safety itself for efficiency
   * reasons. The caller is responsible to ensure this.
   *
   * @param buffer Buffer with a size of 4 complex floats to receive the beam
   * response
   * @param ra Right ascension (rad)
   * @param dec Declination (rad)
   * @param freq Frequency (Hz)
   * @param station_idx Station index
   * @param field_id
   */
  void FullResponse(std::complex<float>* response_matrix, double ra, double dec,
                    double freq, size_t station_idx,
                    size_t field_id) final override;

  aocommon::MC2x2 FullResponse(size_t station_idx, double freq,
                               const vector3r_t& direction,
                               std::mutex* mutex) final override;

  aocommon::MC2x2Diag ArrayFactor(size_t station_idx, double freq,
                                  const vector3r_t& direction,
                                  std::mutex* mutex) final override;

  aocommon::MC2x2 ElementResponse(
      size_t station_idx, double freq,
      const vector3r_t& direction) const final override;

  /**
   * @brief Method for computing the ITRF-vectors, given ra, dec position in
   * radians and using the cached \param time ((MJD(UTC), s))
   */
  void UpdateITRFVectors(double ra, double dec);

 protected:
  casacore::MDirection delay_dir_, tile_beam_dir_;
  vector3r_t station0_, tile0_, dir_itrf_, diff_beam_centre_;
  bool use_channel_frequency_;
  bool use_differential_beam_;
  casacore::MDirection preapplied_beam_dir_;
  CorrectionMode preapplied_correction_mode_;
  double subband_frequency_;

 private:
  /**
   * @brief Update ITRF coordinates for reference station and reference tile
   * direction. Member function leaves the responsibility for providing the
   * mutex to the caller.
   */
  void UpdateITRFVectors(std::mutex& mutex);

  double ra_, dec_;
  std::mutex mutex_;

  // Marks whether the itrf vectors were only partially updated.
  // This bool switches to true if UpdateITRFVectors() is called, since
  // this method doesn't update all the ITRF direction vectors.
  bool has_partial_itrf_update_;
};

}  // namespace pointresponse
}  // namespace everybeam

#endif  // EVERYBEAM_POINTRESPONSE_PHASEDARRAYPOINT_H_
