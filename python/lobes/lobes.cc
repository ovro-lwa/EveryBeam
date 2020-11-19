// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>

#include <Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <iostream>
#include <complex>
#include <H5Cpp.h>
#include <tuple>
#include <string>

#include <boost/math/special_functions/legendre.hpp>

namespace py = pybind11;

// python style enumerate function
// To make it possible to write:
//   for (auto [i, v] : enumerate(iterable)) {...}
template <typename T, typename TIter = decltype(std::begin(std::declval<T>())),
          typename = decltype(std::end(std::declval<T>()))>
constexpr auto enumerate(T &&iterable) {
  struct iterator {
    std::size_t i;
    TIter iter;
    bool operator!=(const iterator &other) const { return iter != other.iter; }
    void operator++() {
      ++i;
      ++iter;
    }
    auto operator*() const { return std::tie(i, *iter); }
  };
  struct iterable_wrapper {
    T iterable;
    auto begin() { return iterator{0, std::begin(iterable)}; }
    auto end() { return iterator{0, std::end(iterable)}; }
  };
  return iterable_wrapper{std::forward<T>(iterable)};
}

double P20(double x) { return 0.5 * (3 * x * x - 1); }
double P21(double x) { return 3.0 * x * std::sqrt(1 - x * x); }
double P22(double x) { return 3 * (1 - x * x); }

int plustwo(int a) {
  std::cout << boost::math::legendre_p(2, 0, 0.5) << '=' << P20(0.5) << '\n'
            << boost::math::legendre_p(2, 1, 0.5) << '=' << P21(0.5) << '\n'
            << boost::math::legendre_p(2, 2, 0.5) << '=' << P22(0.5) << '\n';

  return a + 2;
}

Eigen::MatrixXd big_mat() { return Eigen::MatrixXd::Zero(10000, 10000); }

Eigen::ArrayXd P(int m, int n, py::EigenDRef<Eigen::ArrayXd> x) {
  auto N = x.rows();

  Eigen::ArrayXd result(N);  // = Eigen::VectorXd::Zero(10, 1);
  for (int i = 0; i < N; i++) {
    result[i] = boost::math::legendre_p(n, std::abs(m), x[i]);
  }

  if (m < 0) {
    result *=
        std::pow(-1, -m) * std::tgamma(n + m + 1) / std::tgamma(n - m + 1);
  }

  return std::move(result);
}

Eigen::ArrayXd Pacc(int m, int n, py::EigenDRef<Eigen::ArrayXd> x) {
  // auto N = x.rows();

  return (-(n + m) * (n - m + 1.0) * Eigen::sqrt(1.0 - x * x) * P(m - 1, n, x) -
          m * x * P(m, n, x)) /
         (x * x - 1.0);
}

/** Compute Spherical Wave Function
 *
 * returns std::pair(Etheta, Ephi)
 */
std::pair<Eigen::VectorXcd, Eigen::VectorXcd> F4far_new(
    int s, int m, int n, py::EigenDRef<const Eigen::ArrayXd> theta,
    py::EigenDRef<const Eigen::ArrayXd> phi, double beta) {
  int N = theta.rows();
  Eigen::VectorXd result(N);

  double C;
  if (m) {
    C = beta * std::sqrt(60.0) * 1.0 / std::sqrt(n * (n + 1.0)) *
        std::pow(-m / std::abs(m), m);
  } else {
    C = beta * std::sqrt(60.0) * 1.0 / std::sqrt(n * (n + 1.0));
  }

  Eigen::ArrayXcd q2;
  Eigen::ArrayXcd q3;

  // From cpp >= cpp14, complex literals can be used instead
  constexpr std::complex<double> i_neg = {0.0, -1.0};
  constexpr std::complex<double> i_pos = {0.0, 1.0};

  Eigen::ArrayXd cos_theta = Eigen::cos(theta);
  Eigen::ArrayXd sin_theta = Eigen::sin(theta);
  Eigen::ArrayXd P_cos_theta = P(std::abs(m), n, cos_theta);
  Eigen::ArrayXd Pacc_cos_theta = Pacc(std::abs(m), n, cos_theta);
  Eigen::ArrayXcd exp_i_m_phi = Eigen::exp(i_pos * double(m) * phi);

  if (s == 1) {
    q2 = C * std::pow(i_neg, -n - 1) / beta * i_pos * double(m) /
         (sin_theta)*std::sqrt((2. * n + 1) / 2.0 *
                               std::tgamma(n - std::abs(m) + 1) /
                               std::tgamma(n + std::abs(m) + 1)) *
         P_cos_theta * exp_i_m_phi;

    q3 = C * std::pow(i_neg, -n - 1) / beta *
         std::sqrt((2. * n + 1) / 2.0 * std::tgamma(n - abs(m) + 1) /
                   std::tgamma(n + abs(m) + 1)) *
         Pacc_cos_theta * sin_theta * exp_i_m_phi;
  } else if (s == 2) {
    q2 = -C * std::pow(i_neg, -n) / beta *
         std::sqrt((2. * n + 1) / 2.0 * std::tgamma(n - abs(m) + 1) /
                   std::tgamma(n + abs(m) + 1)) *
         Pacc_cos_theta * sin_theta * exp_i_m_phi;

    q3 = C * std::pow(i_neg, -n) / beta * i_pos * double(m) / sin_theta *
         std::sqrt((2. * n + 1) / 2.0 * std::tgamma(n - abs(m) + 1) /
                   std::tgamma(n + abs(m) + 1)) *
         P_cos_theta * exp_i_m_phi;
  }

  return std::make_pair(q2, q3);
}

PYBIND11_MODULE(pylobes, m) {
  m.doc() = "LOBES module";  // optional module docstring

  m.def("plustwo", &plustwo, "A function which adds two to a number");

  m.def("assoc_legendre",
        static_cast<double (*)(int, double)>(&boost::math::legendre_p),
        "Associated legendre function");

  m.def("scale", [](py::EigenDRef<Eigen::MatrixXd> m, double c) { m *= c; });

  m.def("big_mat", &big_mat, "Big matrix");

  m.def("P", &P, "Legendre");
  m.def("Pacc", &Pacc, "Legendre");
  m.def("F4far_new", &F4far_new, "F4far_new");
}
