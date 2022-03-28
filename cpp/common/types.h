// types.h: Types used in this library.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TYPES_H
#define EVERYBEAM_TYPES_H

#include <array>
#include <cstring>
#include <ostream>
#include <complex>

namespace everybeam {

/** Print the contents of a static array. */
template <typename T, size_t N>
std::ostream& operator<<(std::ostream& out, const std::array<T, N>& obj);

/** Type used for real scalars. */
typedef double real_t;

/** Type used for complex scalars. */
typedef std::complex<double> complex_t;

/** Type used for 2-dimensional real vectors. */
typedef std::array<real_t, 2> vector2r_t;

/** Type used for 3-dimensional real vectors. */
typedef std::array<real_t, 3> vector3r_t;

/** Type used for 2x2 complex diagonal matrices. */
typedef std::array<complex_t, 2> diag22c_t;

/** Type used for 2x2 real matrices. */
typedef std::array<std::array<real_t, 2>, 2> matrix22r_t;

typedef std::array<vector3r_t, 16> TileConfig;

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& out, const std::array<T, N>& obj) {
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
