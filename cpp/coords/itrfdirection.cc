// itrfdirection.cc: Functor that maps time to an ITRF direction.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "itrfdirection.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>

namespace everybeam {
namespace coords {

// ITRF position of CS002LBA, just to use a fixed reference
const vector3r_t ITRFDirection::lofar_position_ = {
    {826577.022720000, 461022.995082000, 5064892.814}};

// TODO: initialize converter with a time (and fixed position) and convert
// specific directions. Needed for wslean as well as for the makeeverybeam
// executable.

ITRFDirection::ITRFDirection(const vector3r_t &position,
                             const vector2r_t &direction) {
  casacore::MVPosition mvPosition(position[0], position[1], position[2]);
  casacore::MPosition mPosition(mvPosition, casacore::MPosition::ITRF);
  frame_ = casacore::MeasFrame(casacore::MEpoch(), mPosition);

  // Order of angles seems to be longitude (along the equator), lattitude
  // (towards the pole).
  casacore::MVDirection mvDirection(direction[0], direction[1]);
  casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);
  converter_ = casacore::MDirection::Convert(
      mDirection,
      casacore::MDirection::Ref(casacore::MDirection::ITRF, frame_));
}

ITRFDirection::ITRFDirection(const vector2r_t &direction)
    : ITRFDirection(lofar_position_, direction) {
  // create ITRF Direction from fixed stationposition
}

ITRFDirection::ITRFDirection(const vector3r_t &position,
                             const vector3r_t &direction) {
  casacore::MVPosition mvPosition(position[0], position[1], position[2]);
  casacore::MPosition mPosition(mvPosition, casacore::MPosition::ITRF);
  frame_ = casacore::MeasFrame(casacore::MEpoch(), mPosition);

  casacore::MVDirection mvDirection(direction[0], direction[1], direction[2]);
  casacore::MDirection mDirection(mvDirection, casacore::MDirection::J2000);
  converter_ = casacore::MDirection::Convert(
      mDirection,
      casacore::MDirection::Ref(casacore::MDirection::ITRF, frame_));
}

ITRFDirection::ITRFDirection(const vector3r_t &direction)
    : ITRFDirection(lofar_position_, direction)

{
  // create ITRF Direction from fixed stationposition
}

vector3r_t ITRFDirection::at(real_t time) const {
  std::lock_guard<std::mutex> lock(itsMutex);

  // Cannot use MeasFrame::resetEpoch(Double), because that assumes the
  // argument is UTC in (fractional) days (MJD).
  frame_.resetEpoch(casacore::Quantity(time, "s"));

  const casacore::MDirection &mITRF = converter_();
  const casacore::MVDirection &mvITRF = mITRF.getValue();

  vector3r_t itrf = {{mvITRF(0), mvITRF(1), mvITRF(2)}};
  return itrf;
}
}  // namespace coords
}  // namespace everybeam