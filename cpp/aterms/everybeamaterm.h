// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_EVERYBEAMATERMS_H_
#define EVERYBEAM_ATERMS_EVERYBEAMATERMS_H_

#include "atermbeam.h"
#include "../beammode.h"

#include <memory>

#include <aocommon/coordinatesystem.h>

namespace casacore {
class MeasurementSet;
}

namespace everybeam {
class Options;

namespace telescope {
class Telescope;
}
namespace aterms {
/**
 * @brief Wraps the EveryBeam::Telescope classes
 * for computing the gridded beam response
 *
 */
class EveryBeamATerm final : public ATermBeam {
 public:
  EveryBeamATerm(const casacore::MeasurementSet& ms,
                 const aocommon::CoordinateSystem& coordinate_system,
                 const everybeam::Options& settings);

 private:
  /**
   * @brief Calculate the gridded response for the \param telescope_
   *
   * @param buffer Buffer
   * @param time Time MJD (s)
   * @param frequency Frequency (Hz)
   * @param fieldId Field id
   */
  bool CalculateBeam(std::complex<float>* buffer, double time, double frequency,
                     size_t fieldId) final override;

  std::unique_ptr<telescope::Telescope> telescope_;
  aocommon::CoordinateSystem coordinate_system_;
  BeamMode beam_mode_;

  size_t cached_field_id;
  double cached_frequency_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
