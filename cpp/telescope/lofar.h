// lofar.h: Base class for computing the response for the LOFAR
// telescope.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_LOFAR_H_
#define EVERYBEAM_TELESCOPE_LOFAR_H_

#include "phasedarray.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <memory>

namespace everybeam {

namespace griddedresponse {
class GriddedResponse;
}  // namespace griddedresponse

namespace pointresponse {
class PointResponse;
}  // namespace pointresponse

namespace telescope {

//! LOFAR telescope class
class [[gnu::visibility("default")]] LOFAR final : public PhasedArray {
 public:
  /**
   * @brief Construct a new LOFAR object
   *
   * @param ms MeasurementSet
   * @param model Element Response model
   * @param options telescope options
   */
  LOFAR(const casacore::MeasurementSet& ms, const Options& options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const aocommon::CoordinateSystem& coordinate_system) const override;

  std::unique_ptr<pointresponse::PointResponse> GetPointResponse(double time)
      const override;

 private:
  bool is_aartfaac_ = false;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_LOFAR_H_
