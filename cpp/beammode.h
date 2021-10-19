// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_BEAMMODE_H_
#define EVERYBEAM_BEAMMODE_H_

#include "correctionmode.h"

namespace everybeam {

/**
 * Describes which beam is computed: array, element or both.
 * These may have different meanings for different telescopes. For LOFAR HBA,
 * the element beam is the dipole beam, and the array beam is the combination
 * of the tile factor and station-array factor.
 *
 * NOTE: \c BeamMode will replace \c CorrectionMode in the future.
 */

using BeamMode = CorrectionMode;

inline BeamMode ParseBeamMode(const std::string& str) {
  return ParseCorrectionMode(str);
}
}  // namespace everybeam

#endif
