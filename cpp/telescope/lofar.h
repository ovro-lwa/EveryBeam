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
class LOFARGrid;
class GriddedResponse;
}  // namespace griddedresponse

namespace pointresponse {
class PointResponse;
class LOFARPoint;
}  // namespace pointresponse

namespace telescope {

//! LOFAR telescope class
class LOFAR final : public PhasedArray {
  friend class griddedresponse::LOFARGrid;
  friend class pointresponse::LOFARPoint;

 public:
  /**
   * @brief Construct a new LOFAR object
   *
   * @param ms MeasurementSet
   * @param model Element Response model
   * @param options telescope options
   */
  LOFAR(const casacore::MeasurementSet &ms, const Options &options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) const override;

  std::unique_ptr<pointresponse::PointResponse> GetPointResponse(
      double time) const override;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_LOFAR_H_
