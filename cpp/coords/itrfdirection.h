// itrfdirection.h: Functor that maps time to an ITRF direction.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ITRFDIRECTION_H
#define EVERYBEAM_ITRFDIRECTION_H

// \file
// Functor that maps time to an ITRF direction.

#include "./../common/types.h"

#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>

#include <memory>
#include <mutex>

namespace everybeam {
namespace coords {
class ITRFDirection {
 public:
  typedef std::shared_ptr<ITRFDirection> Ptr;
  typedef std::shared_ptr<const ITRFDirection> ConstPtr;

  ITRFDirection(const vector3r_t& position, const vector2r_t& direction);
  ITRFDirection(const vector3r_t& position, const vector3r_t& direction);
  ITRFDirection(const vector2r_t& direction);
  ITRFDirection(const vector3r_t& direction);

  vector3r_t at(real_t time) const;

  const static vector3r_t& LOFARPosition() { return lofar_position_; }

 private:
  // ITRF position of CS002LBA, just to use a fixed reference
  const static vector3r_t lofar_position_;

  mutable casacore::MeasFrame frame_;
  mutable casacore::MDirection::Convert converter_;
  mutable std::mutex mutex_;
};
}  // namespace coords
}  // namespace everybeam

#endif
