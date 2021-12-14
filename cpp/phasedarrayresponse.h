// phasedarrayresponse.h: Common functionality for PhasedArrayPoint and
// PhasedArrayGrid.
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_PHASEDARRAYRESPONSE_PHASEDARRAYRESPONSE_H_
#define EVERYBEAM_PHASEDARRAYRESPONSE_PHASEDARRAYRESPONSE_H_

#include "common/types.h"
#include "correctionmode.h"
#include "beammode.h"
#include "beamnormalisationmode.h"

#include <aocommon/matrix2x2.h>
#include <casacore/measures/Measures/MDirection.h>

namespace everybeam {
namespace telescope {

/**
 * @brief Class containing shared functionality and member variables
 * for phased array beamformer responses.
 *
 * TODO: this class might migrate to a better suited namespace in the future.
 */
class PhasedArray;
}  // namespace telescope

class PhasedArrayResponse {
 public:
  ~PhasedArrayResponse() = default;

  /**
   * @brief Determine whether beam normalisation should be applied,
   * and calculate inverse central gain if
   *
   * @param beam_mode Requested beam mode (kNone, kFull, kArrayFactor, kElement)
   * @param time Time (J2000,s)
   * @param frequency Frequency (Hz)
   * @param station_index Station index
   * @param inverse_gain Inverse central gain with which to normalize beam.
   * @return Is beam normalisation required?
   */
  bool CalculateBeamNormalisation(BeamMode beam_mode, double time,
                                  double frequency, size_t station_index,
                                  aocommon::MC2x2F& inverse_gain) const;

 protected:
  PhasedArrayResponse(const telescope::PhasedArray* phasedarray);

  vector3r_t station0_;
  vector3r_t tile0_;
  vector3r_t diff_beam_centre_;

  const casacore::MDirection delay_dir_;
  const casacore::MDirection tile_beam_dir_;
  const casacore::MDirection preapplied_beam_dir_;
  const CorrectionMode preapplied_correction_mode_;

  const BeamNormalisationMode beam_normalisation_mode_;
  // non-const, since OSKAR overwrites these member variables
  // to correct default
  bool use_channel_frequency_;
  double subband_frequency_;

 private:
  const telescope::PhasedArray* phased_array_;
};
}  // namespace everybeam

#endif  // EVERYBEAM_PHASEDARRAYRESPONSE_PHASEDARRAYRESPONSE_H_