// mathutils.h: Various mathematical operations on vectors and matrices.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_MATHUTIL_H
#define EVERYBEAM_MATHUTIL_H

// \file
// Various mathematical operations on vectors and matrices.

#include "types.h"

namespace everybeam {

inline real_t dot(const vector3r_t &arg0, const vector3r_t &arg1) {
  return arg0[0] * arg1[0] + arg0[1] * arg1[1] + arg0[2] * arg1[2];
}

inline double norm(const vector3r_t &arg0) { return sqrt(dot(arg0, arg0)); }

inline vector3r_t operator*(real_t arg0, const vector3r_t arg1) {
  return vector3r_t{arg0 * arg1[0], arg0 * arg1[1], arg0 * arg1[2]};
}

inline vector3r_t operator/(const vector3r_t arg0, real_t arg1) {
  return vector3r_t{arg0[0] / arg1, arg0[1] / arg1, arg0[2] / arg1};
}

inline vector3r_t normalize(const vector3r_t &arg0) {
  return arg0 / norm(arg0);
}

inline vector2r_t cart2thetaphi(const vector3r_t &cart) {
  real_t r = sqrt(cart[0] * cart[0] + cart[1] * cart[1]);
  return vector2r_t{M_PI_2 - atan2(cart[2], r), atan2(cart[1], cart[0])};
}

inline vector3r_t thetaphi2cart(const vector2r_t &thetaphi) {
  real_t r = sin(thetaphi[0]);
  return vector3r_t{r * cos(thetaphi[1]), r * sin(thetaphi[1]),
                    cos(thetaphi[0])};
}

// returns az, el, r.
inline vector3r_t cart2sph(const vector3r_t &cart) {
  real_t r = sqrt(cart[0] * cart[0] + cart[1] * cart[1]);
  return vector3r_t{atan2(cart[1], cart[0]), atan2(cart[2], r), norm(cart)};
}

// expects az, el, r.
inline vector3r_t sph2cart(const vector3r_t &sph) {
  return vector3r_t{sph[2] * cos(sph[1]) * cos(sph[0]),
                    sph[2] * cos(sph[1]) * sin(sph[0]), sph[2] * sin(sph[1])};
}

inline vector3r_t cross(const vector3r_t &arg0, const vector3r_t &arg1) {
  return vector3r_t{arg0[1] * arg1[2] - arg0[2] * arg1[1],
                    arg0[2] * arg1[0] - arg0[0] * arg1[2],
                    arg0[0] * arg1[1] - arg0[1] * arg1[0]};
}

inline vector3r_t operator+(const vector3r_t &arg0, const vector3r_t &arg1) {
  return vector3r_t{arg0[0] + arg1[0], arg0[1] + arg1[1], arg0[2] + arg1[2]};
}

inline vector3r_t operator-(const vector3r_t &arg0, const vector3r_t &arg1) {
  return vector3r_t{arg0[0] - arg1[0], arg0[1] - arg1[1], arg0[2] - arg1[2]};
}
}  // namespace everybeam
#endif  // EVERYBEAM_MATHUTIL_H
