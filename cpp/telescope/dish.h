// dish.h: Base class for dish telescopes (VLA, ATCA, ...).
// Inherits from Telescope class.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_DISH_H_
#define EVERYBEAM_TELESCOPE_DISH_H_

#include "telescope.h"

#include "../circularsymmetric/coefficients.h"

#include <casacore/measures/Measures/MDirection.h>

namespace everybeam {

namespace griddedresponse {
class DishGrid;
}  // namespace griddedresponse

namespace pointresponse {
class PointResponse;
}  // namespace pointresponse

namespace telescope {

/**
 * This class calculates the a-terms for dishes with a circularly symmetric
 * response.
 */
class Dish final : public Telescope {
 public:
  Dish(const casacore::MeasurementSet &ms,
       std::unique_ptr<circularsymmetric::Coefficients> coefficients,
       const Options &options);

  std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) const override;

  std::unique_ptr<pointresponse::PointResponse> GetPointResponse(
      double time) const override;

  /**
   * @brief Get (ra, dec) pointings of fields.
   *
   * @return std::vector<std::pair<double, double>> Vector of size number of
   * fields, and (ra, dec) pointings as entries.
   */
  const std::vector<std::pair<double, double>> &GetFieldPointing() const {
    return field_pointing_;
  }

  const circularsymmetric::Coefficients *GetDishCoefficients() const {
    return dish_coefficients_.get();
  }

 private:
  std::unique_ptr<circularsymmetric::Coefficients> dish_coefficients_;
  /// Store ra, dec pointing per field id from measurement set
  std::vector<std::pair<double, double>> field_pointing_;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_DISH_H_
