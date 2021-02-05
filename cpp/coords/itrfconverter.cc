// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Dir2ITRF.cc: Convertor that maps time to an ITRF direction.

#include "itrfdirection.h"
#include "itrfconverter.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>

namespace everybeam {
namespace coords {
// TODO: Initialize converter with a time (and fixed position) and convert
// specific directions.
//      Needed for wslean as well as for the makeeverybeam executable.

ITRFConverter::ITRFConverter(real_t time) {
  // create ITRF Direction from fixed stationposition
  casacore::MVPosition mv_position(ITRFDirection::LOFARPosition()[0],
                                   ITRFDirection::LOFARPosition()[1],
                                   ITRFDirection::LOFARPosition()[2]);
  casacore::MPosition m_position(mv_position, casacore::MPosition::ITRF);
  casacore::MEpoch time_epoch(casacore::Quantity(time, "s"));
  frame_ = casacore::MeasFrame(time_epoch, m_position);

  // Order of angles seems to be longitude (along the equator), lattitude
  // (towards the pole).
  converter_ = casacore::MDirection::Convert(
      casacore::MDirection::J2000,
      casacore::MDirection::Ref(casacore::MDirection::ITRF, frame_));
}

void ITRFConverter::SetTime(real_t time) {
  // Cannot use MeasFrame::resetEpoch(Double), because that assumes the
  // argument is UTC in (fractional) days (MJD).
  frame_.resetEpoch(casacore::Quantity(time, "s"));
}

vector3r_t ITRFConverter::j2000ToITRF(const vector2r_t &j2000Direction) const {
  casacore::MVDirection mv_direction(j2000Direction[0], j2000Direction[1]);
  casacore::MDirection m_direction(mv_direction, casacore::MDirection::J2000);
  const casacore::MVDirection mvITRF = converter_(m_direction).getValue();

  return vector3r_t{{mvITRF(0), mvITRF(1), mvITRF(2)}};
}

vector3r_t ITRFConverter::j2000ToITRF(const vector3r_t &j2000Direction) const {
  casacore::MVDirection mv_direction(j2000Direction[0], j2000Direction[1],
                                     j2000Direction[2]);
  casacore::MDirection m_direction(mv_direction, casacore::MDirection::J2000);

  const casacore::MVDirection mvITRF = converter_(m_direction).getValue();

  return vector3r_t{{mvITRF(0), mvITRF(1), mvITRF(2)}};
}

vector3r_t ITRFConverter::ToITRF(const casacore::MDirection &direction) const {
  const casacore::MVDirection mvITRF = converter_(direction).getValue();
  return vector3r_t{{mvITRF(0), mvITRF(1), mvITRF(2)}};
}

casacore::MDirection ITRFConverter::ToDirection(
    const vector2r_t &j2000Direction) const {
  casacore::MVDirection mv_direction(j2000Direction[0], j2000Direction[1]);
  casacore::MDirection m_direction(mv_direction, casacore::MDirection::J2000);
  return converter_(m_direction);
}

casacore::MDirection ITRFConverter::ToDirection(
    const vector3r_t &j2000Direction) const {
  casacore::MVDirection mv_direction(j2000Direction[0], j2000Direction[1],
                                     j2000Direction[2]);
  casacore::MDirection m_direction(mv_direction, casacore::MDirection::J2000);

  return converter_(m_direction);
}

casacore::MDirection ITRFConverter::ToDirection(
    const casacore::MDirection &direction) const {
  return converter_(direction);
}
}  // namespace coords
}  // namespace everybeam