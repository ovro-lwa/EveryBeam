// mwa.h: Class for MWA telescopes.
// Inherits from Telescope class.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_MWA_H_
#define EVERYBEAM_TELESCOPE_MWA_H_

#include "telescope.h"

#include <casacore/measures/Measures/MPosition.h>

namespace everybeam {

namespace griddedresponse {
class GriddedResponse;
}  // namespace griddedresponse

namespace pointresponse {
class PointResponse;
}  // namespace pointresponse

namespace telescope {

class MWA final : public Telescope {
 public:
  /**
   * @brief Construct a new MWA object
   *
   * @param ms MeasurementSet
   * @param model Element Response model
   * @param options telescope options
   */
  MWA(const casacore::MeasurementSet &ms, const Options &options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) const override;

  std::unique_ptr<pointresponse::PointResponse> GetPointResponse(
      double time) const override;

  casacore::MPosition GetArrayPosition() const { return array_position_; }
  const std::array<double, 16> &GetDelays() const { return delays_; }

 private:
  casacore::MPosition array_position_;
  std::array<double, 16> delays_;
};

}  // namespace telescope
}  // namespace everybeam
#endif  // EVERYBEAM_TELESCOPE_MWA_H_
