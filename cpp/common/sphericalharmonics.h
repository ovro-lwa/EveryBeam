// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_COMMON_SPHERICAL_HARMONICS_H_
#define EVERYBEAM_COMMON_SPHERICAL_HARMONICS_H_

#include <Eigen/Core>

namespace everybeam {
namespace common {

/**
 * @brief Evaluate associated Legendre polynomial
 *
 * @param m the order of the polynomial
 * @param n the degree of the polynomial
 * @param x the argument
 * @return double
 */
[[gnu::visibility("default")]] double P(int m, int n, double x);

/**
 * @brief Evaluate associated Legendre polynomial (array input)
 */
[[gnu::visibility("default")]] Eigen::ArrayXd P(int m, int n,
                                                const Eigen::ArrayXd& x);

/**
 * @brief P' of associated Legendre polynomial
 *
 * @param m the order of the polynomial
 * @param n the degree of the polynomial
 * @param x the argument
 * @return double
 */
[[gnu::visibility("default")]] double Pacc(int m, int n, double x);

/**
 * @brief P' of asscociated Legendre polynomial (array input)
 */
[[gnu::visibility("default")]] Eigen::ArrayXd Pacc(int m, int n,
                                                   const Eigen::ArrayXd& x);

/**
 * @brief Evaluate spherical harmonic for given position
 *
 * @param s
 * @param m the order of the Legendre polynomial
 * @param n the degree of the Legendre polynomial
 * @param theta zenith angle [rad]
 * @param phi elevation angle [rad]
 * @return std::pair<std::complex<double>, std::complex<double>>
 */
[[gnu::visibility(
    "default")]] std::pair<std::complex<double>, std::complex<double>>
F4far_new(int s, int m, int n, double theta, double phi);

/**
 * @brief Evaluate spherical harmonics for given positions (vector input)
 */
[[gnu::visibility("default")]] std::pair<Eigen::VectorXcd, Eigen::VectorXcd>
F4far_new(int s, int m, int n, const Eigen::VectorXd& theta,
          const Eigen::VectorXd& phi);

}  // namespace common
}  // namespace everybeam
#endif