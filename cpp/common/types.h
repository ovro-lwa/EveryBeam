// types.h: Types used in this library.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TYPES_H
#define EVERYBEAM_TYPES_H

// \file
// Types used in this library.

#include <array>
#include <cstring>
#include <ostream>
#include <complex>

namespace everybeam {

/** Print the contents of a static array. */
template <typename T, size_t N>
std::ostream &operator<<(std::ostream &out, const std::array<T, N> &obj);

/** Type used for real scalars. */
typedef double real_t;

/** Type used for complex scalars. */
typedef std::complex<double> complex_t;

/** Type used for 2-dimensional real vectors. */
typedef std::array<real_t, 2> vector2r_t;

/** Type used for 3-dimensional real vectors. */
typedef std::array<real_t, 3> vector3r_t;

/** Type used for 2x2 real diagonal matrices. */
typedef std::array<real_t, 2> diag22r_t;

/** Type used for 2x2 complex diagonal matrices. */
typedef std::array<complex_t, 2> diag22c_t;

/** Type used for 2x2 real matrices. */
typedef std::array<std::array<real_t, 2>, 2> matrix22r_t;

/** Type used for 2x2 complex matrices. */
typedef std::array<std::array<complex_t, 2>, 2> matrix22c_t;

/** Response of an array of antenna elements. */
struct raw_response_t {
  /** Combined response of all (enabled) antenna elements in the array. */
  matrix22c_t response;

  /** Number of antenna elements contributing to the combined response, per
   *  polarization.
   */
  diag22r_t weight;
};

/** Array factor of an array of antenna elements. A wave front of an incident
 *  plane wave will arrive at each antenna element at a potentially different
 *  time. The time of arrival depends on the location of the antenna element and
 *  the direction of arrival of the plane wave. With respect to a pre-defined
 *  phase reference location, there is a (possibly negative) delay between the
 *  arrival of a wave front at a given antenna element and the arrival of the
 *  same wave front at the phase reference location. The array factor is the sum
 *  of the phase shifts due to these delays. It describes the "sensitivity" of
 *  the array as a function of direction.
 */
struct raw_array_factor_t {
  /** Array factor due to all (enabled) antenna elements in the array. */
  diag22c_t factor;

  /** Number of antenna elements contributing to the array factor, per
   *  polarization.
   */
  diag22r_t weight;
};

typedef std::array<vector3r_t, 16> TileConfig;

template <typename T, size_t N>
std::ostream &operator<<(std::ostream &out, const std::array<T, N> &obj) {
  // print(out, obj.begin(), obj.end());
  out << "[";
  for (auto it : obj) {
    out << it;
    if (it != *obj.rbegin()) out << ", ";
  }
  out << "]\n";
  return out;
}

}  // namespace everybeam

#endif  // EVERYBEAM_TYPES_H
