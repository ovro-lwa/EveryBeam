// aartfaacgrid.h: Class for computing the gridded response for
// LOFAR observations in AARTFAAC modus
//
// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_AARTFAACGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_AARTFAACGRID_H_

#include "phasedarraygrid.h"
#include "../telescope/telescope.h"

namespace everybeam {
namespace griddedresponse {

/**
 * @brief Class for computing the AARTFAAC beam responses on a
 * user specified grid. Class exploits the
 * fact that the station response is identical for all stations
 * within an AARTFAAC observation.
 */
class [[gnu::visibility("default")]] AartfaacGrid final
    : public PhasedArrayGrid {
 public:
  /**
   * @brief Construct a new AartfaacGrid object
   *
   * @param telescope_ptr Pointer to telescope::LOFAR object
   * @param coordinate_system CoordinateSystem struct
   */
  AartfaacGrid(const telescope::Telescope* telescope_ptr,
               const aocommon::CoordinateSystem& coordinate_system)
      : PhasedArrayGrid(telescope_ptr, coordinate_system){};

  void ResponseAllStations(BeamMode beam_mode, std::complex<float> * buffer,
                           double time, double frequency, size_t field_id)
      override;

 private:
  // Override MakeIntegratedSnapshot for efficiency.
  // Implementation is similar to MWAGrid::MakeIntegratedSnapshot
  void MakeIntegratedSnapshot(BeamMode beam_mode,
                              std::vector<aocommon::HMC4x4> & matrices,
                              double time, double frequency, size_t field_id,
                              const double* baseline_weights_interval) override;
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_AARTFAACGRID_H_
