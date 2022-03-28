// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_SKAMID_H_
#define EVERYBEAM_TELESCOPE_SKAMID_H_

#include "telescope.h"

namespace everybeam {

namespace griddedresponse {
class GriddedResponse;
}

namespace pointresponse {
class PointResponse;
}

namespace telescope {
//! SKA-MID telescope class
class SkaMid final : public Telescope {
 public:
  SkaMid(const casacore::MeasurementSet& ms, const Options& options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem& coordinate_system) const override;

  std::unique_ptr<pointresponse::PointResponse> GetPointResponse(
      double time) const override;

  /**
   * @brief Get (ra, dec) pointings of fields.
   *
   * @return std::vector<std::pair<double, double>> Vector of size number of
   * fields, and (ra, dec) pointings as entries.
   */
  const std::vector<std::pair<double, double>>& GetFieldPointing() const {
    return field_pointing_;
  }

  /**
   * @brief Diameter of SKA-MID dish (m)
   */
  double GetDiameter() const;

  /**
   * @brief Blockage of SKA-MID due to receiver (m)
   */
  double GetBlockage() const;

 private:
  std::vector<std::pair<double, double>> field_pointing_;

  ElementResponseModel element_response_model_;
};
}  // namespace telescope
}  // namespace everybeam
#endif