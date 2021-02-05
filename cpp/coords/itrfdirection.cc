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
  casacore::MVPosition mv_position(position[0], position[1], position[2]);
  casacore::MPosition m_position(mv_position, casacore::MPosition::ITRF);
  frame_ = casacore::MeasFrame(casacore::MEpoch(), m_position);

  // Order of angles seems to be longitude (along the equator), lattitude
  // (towards the pole).
  casacore::MVDirection mv_direction(direction[0], direction[1]);
  casacore::MDirection m_direction(mv_direction, casacore::MDirection::J2000);
  converter_ = casacore::MDirection::Convert(
      m_direction,
      casacore::MDirection::Ref(casacore::MDirection::ITRF, frame_));
}

ITRFDirection::ITRFDirection(const vector2r_t &direction)
    : ITRFDirection(lofar_position_, direction) {
  // create ITRF Direction from fixed stationposition
}

ITRFDirection::ITRFDirection(const vector3r_t &position,
                             const vector3r_t &direction) {
  casacore::MVPosition mv_position(position[0], position[1], position[2]);
  casacore::MPosition m_position(mv_position, casacore::MPosition::ITRF);
  frame_ = casacore::MeasFrame(casacore::MEpoch(), m_position);

  casacore::MVDirection mv_direction(direction[0], direction[1], direction[2]);
  casacore::MDirection m_direction(mv_direction, casacore::MDirection::J2000);
  converter_ = casacore::MDirection::Convert(
      m_direction,
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

  const casacore::MDirection &m_ITRF = converter_();
  const casacore::MVDirection &mv_ITRF = m_ITRF.getValue();

  vector3r_t itrf = {{mv_ITRF(0), mv_ITRF(1), mv_ITRF(2)}};
  return itrf;
}
}  // namespace coords
}  // namespace everybeam