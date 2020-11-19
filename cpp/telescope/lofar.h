// lofar.h: Base class for computing the response for the LOFAR
// telescope.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_LOFAR_H_
#define EVERYBEAM_TELESCOPE_LOFAR_H_

#include "../station.h"
#include "../elementresponse.h"
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

namespace telescope {

//! LOFAR telescope class
class LOFAR final : public PhasedArray {
  friend class griddedresponse::LOFARGrid;

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
      const coords::CoordinateSystem &coordinate_system) override;

 private:
  Station::Ptr ReadStation(const casacore::MeasurementSet &ms,
                           const std::size_t id,
                           const ElementResponseModel model) const override;

  std::vector<Station::Ptr> stations_;
};
}  // namespace telescope
}  // namespace everybeam

#endif  // EVERYBEAM_TELESCOPE_LOFAR_H_
