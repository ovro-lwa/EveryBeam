// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "klfittingaterm.h"
#include <aocommon/imagecoordinates.h>

namespace everybeam {
namespace aterms {

KlFittingATerm::KlFittingATerm(
    const coords::CoordinateSystem& coordinate_system)
    : coordinate_system_(coordinate_system) {}

bool KlFittingATerm::Calculate(std::complex<float>* buffer, double time,
                               double frequency, size_t field_id,
                               const double* uvw_in_m) {
  throw std::runtime_error("KlFittingATerm::Calculate not implemented.");
  return false;
}

double KlFittingATerm::AverageUpdateTime() const {
  throw std::runtime_error(
      "KlFittingATerm::AverageUpdateTime not implemented.");
  return 0.0;
}

}  // namespace aterms
}  // namespace everybeam
