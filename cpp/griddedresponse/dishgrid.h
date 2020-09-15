// dishgrid.h: Class for computing the circular symmetric (gridded) response.
//
// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the EveryBeam software suite.
// The EveryBeam software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The EveryBeam software suite is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the EveryBeam software suite. If not, see
// <http://www.gnu.org/licenses/>.
//
// $Id$

#ifndef EVERYBEAM_GRIDDEDRESPONSE_DISHGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_DISHGRID_H_

#include "griddedresponse.h"

namespace everybeam {
namespace griddedresponse {

/**
 * @brief Class for computing the gridded response of dish telescopes,
 * e.g. VLA, ATCA.
 *
 */
class DishGrid final : public GriddedResponse {
 public:
  DishGrid(telescope::Telescope* telescope_ptr,
           const coords::CoordinateSystem coordinate_system)
      : GriddedResponse(telescope_ptr, coordinate_system){};

  void CalculateStation(std::complex<float>* buffer, double time,
                        double frequency, size_t station_idx,
                        size_t field_id) override;

  void CalculateAllStations(std::complex<float>* buffer, double time,
                            double frequency, size_t field_id) override;

  virtual void CalculateIntegratedResponse(
      double* buffer, double time, double frequency, size_t field_id,
      size_t undersampling_factor,
      const std::vector<double>& baseline_weights) override;

  virtual void CalculateIntegratedResponse(
      double* buffer, const std::vector<double>& time_array, double frequency,
      size_t field_id, size_t undersampling_factor,
      const std::vector<double>& baseline_weights) override {
    // Time does not play a role in the integrated response of a dish telescope,
    // so call CalculateIntegratedResponse as if it were one time step
    CalculateIntegratedResponse(buffer, 0., frequency, field_id,
                                undersampling_factor, std::vector<double>{0});
  };

 private:
  /**
   * @brief Make integrated snapshot, specialized/simplified for dish
   * telescopes.
   *
   * @param matrices Vector of Mueller matrices
   * @param frequency Frequency (Hz)
   * @param field_id Field id
   */
  void MakeIntegratedSnapshot(std::vector<aocommon::HMC4x4>& matrices,
                              double frequency, size_t field_id);
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_DISHGRID_H_