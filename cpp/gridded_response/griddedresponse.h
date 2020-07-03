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

#include <memory>

namespace everybeam {

namespace telescope {
class Telescope;
}

// TODO: just temporary location!
struct CoordinateSystem {
  std::size_t width, height;
  double ra, dec, dl, dm, phaseCentreDL, phaseCentreDM;
};

namespace gridded_response {

/**
 * @brief Virtual base class to compute the gridded response
 *
 */
class GriddedResponse {
 public:
  typedef std::unique_ptr<GriddedResponse> Ptr;

  // TODO: can be deprecated in a later stage
  virtual void CalculateStation(std::size_t station_id) = 0;
  virtual void CalculateStation() = 0;
  virtual void CalculateAllStations() = 0;

  // TODO: complete!
  virtual void CalculateStation(std::complex<float>* buffer, size_t station_id,
                                double time, double freq) = 0;
  // Repeated call of calculate single?
  virtual void CalculateAllStations(std::complex<float>* buffer, double time,
                                    double freq) = 0;

 protected:
  /**
   * @brief Construct a new Gridded Response object
   *
   * @param telescope_ptr Pointer to telescope::Telescope object
   * @param coordinateSystem CoordinateSystem struct
   */
  GriddedResponse(const telescope::Telescope* telescope_ptr,
                  const CoordinateSystem& coordinateSystem)
      : _telescope(telescope_ptr),
        _width(coordinateSystem.width),
        _height(coordinateSystem.height),
        _ra(coordinateSystem.ra),
        _dec(coordinateSystem.dec),
        _dl(coordinateSystem.dl),
        _dm(coordinateSystem.dm),
        _phaseCentreDL(coordinateSystem.phaseCentreDL),
        _phaseCentreDM(coordinateSystem.phaseCentreDM){};

  const telescope::Telescope* _telescope;
  size_t _width, _height;
  double _ra, _dec, _dl, _dm, _phaseCentreDL, _phaseCentreDM;
};
}  // namespace gridded_response
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_GRIDDEDRESPONSE_H_