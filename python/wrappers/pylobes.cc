// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "../../cpp/common/sphericalharmonics.h"

namespace py = pybind11;

using everybeam::common::F4far_new;
using everybeam::common::P;
using everybeam::common::Pacc;

void init_lobes(py::module& m) {
  py::module m_sub = m.def_submodule("lobes");

  m_sub.doc() = "Helper methods for LOBEs response model";

  m_sub.def("P", static_cast<double (*)(int, int, double)>(&P),
            "Evaluate Legendre polynomial (scalar input)");
  m_sub.def(
      "P",
      [](int m, int n, py::EigenDRef<const Eigen::ArrayXd> x) {
        return P(m, n, x);
      },
      "Evaluate Legendre polynomial (array input)");

  // TODO: correct naming for Pacc?!
  m_sub.def("Pacc", static_cast<double (*)(int, int, double)>(&Pacc),
            "Evaluate Legendre polynomial (scalar input)");
  m_sub.def(
      "Pacc",
      [](int m, int n, py::EigenDRef<const Eigen::ArrayXd> x) {
        return Pacc(m, n, x);
      },
      "Evaluate Legendre polynomial (array input)");

  m_sub.def(
      "F4far_new",
      static_cast<std::pair<std::complex<double>, std::complex<double>> (*)(
          int, int, int, double, double)>(&F4far_new),
      "Evaluate spherical harmonics at given (theta, phi) angle");
  m_sub.def("F4far_new",
            static_cast<std::pair<Eigen::VectorXcd, Eigen::VectorXcd> (*)(
                int, int, int, const Eigen::VectorXd&, const Eigen::VectorXd&)>(
                &everybeam::common::F4far_new),
            "Evaluate spherical harmonics at given theta/phi angles "
            "(vector/flattened array input)");
}