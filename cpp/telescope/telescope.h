// telescope.h: Base class for computing the Telescope response.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_TELESCOPE_H_
#define EVERYBEAM_TELESCOPE_TELESCOPE_H_

#include "../options.h"

#include <vector>
#include <memory>
#include <cassert>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>

namespace everybeam {

namespace griddedresponse {
class GriddedResponse;
}

namespace coords {
struct CoordinateSystem;
}

namespace telescope {

/**
 * @brief Telescope class, forming the base class for specific telescopes.
 *
 */
class Telescope {
 public:
  virtual ~Telescope(){};

  /**
   * @brief Return the gridded response object
   *
   * @param coordinate_system Coordinate system struct
   * @return GriddedResponse::Ptr
   */
  virtual std::unique_ptr<griddedresponse::GriddedResponse> GetGriddedResponse(
      const coords::CoordinateSystem &coordinate_system) = 0;

  bool GetIsTimeRelevant() const { return is_time_relevant_; };
  std::size_t GetNrStations() const { return nstations_; };
  Options GetOptions() const { return options_; };

 protected:
  /**
   * @brief Construct a new Telescope object
   *
   * @param ms MeasurementSet
   * @param options telescope options
   */
  Telescope(const casacore::MeasurementSet &ms, const Options &options)
      : nstations_(ms.antenna().nrow()), options_(options){};

  void SetIsTimeRelevant(bool is_time_relevant) {
    is_time_relevant_ = is_time_relevant;
  };

  std::size_t nstations_;
  Options options_;

 private:
  bool is_time_relevant_ = true;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_TELESCOPE_H_
