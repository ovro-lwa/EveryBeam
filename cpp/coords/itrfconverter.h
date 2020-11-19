// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// itrfconverter.h: Convert time to an ITRF direction.

#ifndef EVERYBEAM_DIR2ITRF_H
#define EVERYBEAM_DIR2ITRF_H

// \file
// Functor that maps J2000 to an ITRF direction.

#include "./../common/types.h"

#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>

#include <memory>

namespace everybeam {
namespace coords {
/**
 * @brief Class providing utilities for coordinate transformations
 * to and from ITRF (International Terrestrial Reference Frame).
 * NOTE: this class is not thread-safe due to casacore dependencies.
 *
 */
class ITRFConverter {
 public:
  ITRFConverter(real_t time);

  void SetTime(real_t time);
  vector3r_t j2000ToITRF(const vector2r_t &j2000Direction) const;
  vector3r_t j2000ToITRF(const vector3r_t &j2000Direction) const;
  vector3r_t ToITRF(const casacore::MDirection &direction) const;
  casacore::MDirection ToDirection(const vector2r_t &j2000Direction) const;
  casacore::MDirection ToDirection(const vector3r_t &j2000Direction) const;
  casacore::MDirection ToDirection(const casacore::MDirection &direction) const;

 private:
  casacore::MeasFrame frame_;
  mutable casacore::MDirection::Convert converter_;
};
}  // namespace coords
}  // namespace everybeam
#endif
