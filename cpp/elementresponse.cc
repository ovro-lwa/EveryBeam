// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "elementresponse.h"

namespace everybeam {
std::ostream& operator<<(std::ostream& os, ElementResponseModel model) {
  switch (model) {
    case kDefault:
      os << "Default";
      break;
    case kHamaker:
      os << "Hamaker";
      break;
    case kLOBES:
      os << "LOBES";
      break;
    case kOSKARDipole:
      os << "OSKARDipole";
      break;
    case kOSKARSphericalWave:
      os << "OSKARSphericalWave";
      break;
    default:
      os.setstate(std::ios_base::failbit);
  }
  return os;
}
}  // namespace everybeam
