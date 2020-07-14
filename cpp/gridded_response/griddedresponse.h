// GriddedResponse.h: Base class for computing the (gridded) telescope response.
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

#ifndef EVERYBEAM_GRIDDEDRESPONSE_GRIDDEDRESPONSE_H_
#define EVERYBEAM_GRIDDEDRESPONSE_GRIDDEDRESPONSE_H_

#include "./../coords/coord_utils.h"
#include "./../coords/ITRFDirection.h"
#include "./../coords/ITRFConverter.h"

#include <memory>
#include <vector>
#include <thread>
#include <aocommon/lane.h>
#include <casacore/measures/Measures/MDirection.h>

namespace everybeam {

namespace telescope {
class Telescope;
}

namespace gridded_response {

/**
 * @brief Virtual base class to compute the gridded response
 *
 */
class GriddedResponse {
 public:
  /**
   * @brief Compute the Beam response for a single station
   *
   * @param buffer Output buffer, compute and set size with
   * GriddedResponse::GetBufferSize(1)
   * @param station_idx Station index, must be smaller than number of stations
   * in the Telescope
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param frequency Frequency (Hz)
   */
  virtual bool CalculateStation(std::complex<float>* buffer, double time,
                                double freq, const size_t station_idx) = 0;

  /**
   * @brief Compute the Beam response for all stations in a Telescope
   *
   * @param buffer Output buffer, compute and set size with
   * GriddedResponse::GetBufferSize()
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param frequency Frequency (Hz)
   */
  virtual bool CalculateAllStations(std::complex<float>* buffer, double time,
                                    double frequency) = 0;

  std::size_t GetBufferSize(std::size_t nstations) {
    return std::size_t(nstations * width_ * height_ * 2 * 2);
  }

 protected:
  /**
   * @brief Construct a new Gridded Response object
   *
   * @param telescope_ptr Pointer to telescope::Telescope object
   * @param coordinate_system CoordinateSystem struct
   */
  GriddedResponse(const telescope::Telescope* const telescope_ptr,
                  const coords::CoordinateSystem& coordinate_system)
      : telescope_(telescope_ptr),
        width_(coordinate_system.width),
        height_(coordinate_system.height),
        ra_(coordinate_system.ra),
        dec_(coordinate_system.dec),
        dl_(coordinate_system.dl),
        dm_(coordinate_system.dm),
        phase_centre_dl_(coordinate_system.phase_centre_dl),
        phase_centre_dm_(coordinate_system.phase_centre_dm){};

  const telescope::Telescope* const telescope_;
  size_t width_, height_;
  double ra_, dec_, dl_, dm_, phase_centre_dl_, phase_centre_dm_;
};
}  // namespace gridded_response
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_GRIDDEDRESPONSE_H_