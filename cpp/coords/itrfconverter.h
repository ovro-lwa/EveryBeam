// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
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

namespace everybeam {
namespace coords {
/**
 * @brief Class providing utilities for coordinate transformations
 * to and from ITRF (International Terrestrial Reference Frame).
 * NOTE: this class is not thread-safe due to casacore dependencies.
 *
 */
class ItrfConverter {
 public:
  ItrfConverter(real_t time);

  vector3r_t RaDecToItrf(double ra, double dec) const;
  vector3r_t ToItrf(const casacore::MDirection& direction) const;

 private:
  casacore::MeasFrame frame_;
  mutable casacore::MDirection::Convert converter_;
};
}  // namespace coords
}  // namespace everybeam
#endif
