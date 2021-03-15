// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "sphericalharmonics.h"

#include <cmath>
#include <Eigen/Core>
#include <boost/math/special_functions/legendre.hpp>

namespace everybeam {
namespace common {

namespace {
template <typename T>
T PaccTemplate(int m, int n, const T& x, const T& sqrt_x) {
  return (-(n + m) * (n - m + 1.0) * sqrt_x * P(m - 1, n, x) -
          m * x * P(m, n, x)) /
         (x * x - 1.0);
}

template <typename T, typename U>
void F4farNewTemplate(int s, int m, int n, const T& cos_theta,
                      const T& sin_theta, const U& exp_phi, U& q2, U& q3) {
  double C = std::sqrt(60.0) / std::sqrt(n * (n + 1.0));
  if (m) {
    C *= std::pow(-m / std::abs(m), m);
  }

  // From cpp >= cpp14, complex literals can be used instead
  constexpr std::complex<double> i_neg = {0.0, -1.0};
  constexpr std::complex<double> i_pos = {0.0, 1.0};

  T P_cos_theta = P(std::abs(m), n, cos_theta);
  T Pacc_cos_theta = Pacc(std::abs(m), n, cos_theta);

  if (s == 1) {
    q2 = C * std::pow(i_neg, -n - 1) * i_pos * double(m) /
         (sin_theta)*std::sqrt((2. * n + 1) / 2.0 *
                               std::tgamma(n - std::abs(m) + 1) /
                               std::tgamma(n + std::abs(m) + 1)) *
         P_cos_theta * exp_phi;

    q3 = C * std::pow(i_neg, -n - 1) *
         std::sqrt((2. * n + 1) / 2.0 * std::tgamma(n - abs(m) + 1) /
                   std::tgamma(n + abs(m) + 1)) *
         Pacc_cos_theta * sin_theta * exp_phi;
  } else if (s == 2) {
    q2 = -C * std::pow(i_neg, -n) *
         std::sqrt((2. * n + 1) / 2.0 * std::tgamma(n - abs(m) + 1) /
                   std::tgamma(n + abs(m) + 1)) *
         Pacc_cos_theta * sin_theta * exp_phi;

    q3 = C * std::pow(i_neg, -n) * i_pos * double(m) / sin_theta *
         std::sqrt((2. * n + 1) / 2.0 * std::tgamma(n - abs(m) + 1) /
                   std::tgamma(n + abs(m) + 1)) *
         P_cos_theta * exp_phi;
  }
}
}  // namespace

double P(int m, int n, double x) {
  double result = boost::math::legendre_p(n, std::abs(m), x);
  if (m < 0) {
    int phase = ((-m) % 2 == 0) ? 1 : -1;
    result *= phase * std::tgamma(n + m + 1) / std::tgamma(n - m + 1);
  }
  return result;
}

Eigen::ArrayXd P(int m, int n, const Eigen::ArrayXd& x) {
  const int N = x.rows();
  Eigen::ArrayXd result(N);
  for (int i = 0; i < N; i++) {
    result[i] = boost::math::legendre_p(n, std::abs(m), x[i]);
  }

  if (m < 0) {
    int phase = ((-m) % 2 == 0) ? 1 : -1;
    result *= phase * std::tgamma(n + m + 1) / std::tgamma(n - m + 1);
  }
  return result;
}

double Pacc(int m, int n, double x) {
  return PaccTemplate<double>(m, n, x, std::sqrt(1.0 - x * x));
}

Eigen::ArrayXd Pacc(int m, int n, const Eigen::ArrayXd& x) {
  return PaccTemplate<Eigen::ArrayXd>(m, n, x, Eigen::sqrt(1.0 - x * x));
}

std::pair<std::complex<double>, std::complex<double>> F4far_new(int s, int m,
                                                                int n,
                                                                double theta,
                                                                double phi) {
  std::complex<double> q2, q3;
  // Avoid singularities due to theta = 0
  if (std::abs(theta) < 1e-6) {
    theta = 1e-6;
  }
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  std::complex<double> exp_i_m_phi =
      exp(std::complex<double>{0.0, 1.0} * double(m) * phi);

  F4farNewTemplate(s, m, n, cos_theta, sin_theta, exp_i_m_phi, q2, q3);
  return std::make_pair(q2, q3);
}

std::pair<Eigen::VectorXcd, Eigen::VectorXcd> F4far_new(
    int s, int m, int n, const Eigen::VectorXd& theta,
    const Eigen::VectorXd& phi) {
  Eigen::ArrayXcd q2;
  Eigen::ArrayXcd q3;

  // Avoid singularities due to theta = 0
  auto theta_in = theta;
  theta_in =
      theta.unaryExpr([](double v) { return std::abs(v) >= 1e-6 ? v : 1e-6; });

  const Eigen::ArrayXd cos_theta = Eigen::cos(theta_in.array());
  const Eigen::ArrayXd sin_theta = Eigen::sin(theta_in.array());
  const Eigen::ArrayXcd exp_i_m_phi =
      Eigen::exp(std::complex<double>{0.0, 1.0} * double(m) * phi.array());

  F4farNewTemplate(s, m, n, cos_theta, sin_theta, exp_i_m_phi, q2, q3);
  return std::make_pair(q2, q3);
}

}  // namespace common
}  // namespace everybeam