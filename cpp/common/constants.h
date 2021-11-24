// constants.h: (numeric/string) constants used in this library.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_CONSTANTS_H
#define EVERYBEAM_CONSTANTS_H

#include "types.h"

namespace everybeam {
namespace common {
/** Speed of light (m/s) */
constexpr real_t c = 2.99792458e+08;

/** 0.25*pi */
constexpr real_t pi_4 = 0.7853981633974483096156608;
}  // namespace common
}  // namespace everybeam

#endif
