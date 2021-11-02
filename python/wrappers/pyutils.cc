// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "elementresponse.h"
#include "options.h"
#include "common/mathutils.h"
#include "coords/coordutils.h"

namespace py = pybind11;
using everybeam::cart2thetaphi;
using everybeam::ElementResponseModel;
using everybeam::Options;
using everybeam::thetaphi2cart;
using everybeam::vector2r_t;
using everybeam::vector3r_t;
using everybeam::coords::CoordinateSystem;

namespace {
// Convert pyarray of size 3 to vector3r_t
vector3r_t np2vector3r_t(const py::array_t<double> pyarray) {
  auto r = pyarray.unchecked<1>();
  if (r.size() != 3) {
    throw std::runtime_error("Pyarry is of incorrect size, must be 3.");
  }
  return vector3r_t{r[0], r[1], r[2]};
}
}  // namespace

void init_utils(py::module &m) {
  m.def(
      "cart2thetaphi",
      [](py::array_t<double> pydirection) {
        vector3r_t direction = np2vector3r_t(pydirection);
        vector2r_t thetaphi = cart2thetaphi(direction);
        return thetaphi;
      },
      R"pbdoc(
        Convert direction vector (ITRF or local East-North-Up)
        to theta, phi.

        Parameters
        ----------
        direction: np.1darray
            Direction vector

        Returns
        -------
        list:
            [theta, phi]
       )pbdoc",
      py::arg("direction"));
  m.def(
      "thetaphi2cart",
      [](double theta, double phi) -> py::array_t<float> {
        vector2r_t thetaphi = {theta, phi};
        vector3r_t direction = everybeam::thetaphi2cart(thetaphi);
        return py::array_t<float>(py::cast(direction));
      },
      R"pbdoc(
        Convert theta, phi angles to direction vector

        Parameters
        ----------
        theta: float
            theta angle [rad]
        phi: float
            phi angle [rad]

        Returns
        -------
        np.1darray:

       )pbdoc",
      py::arg("theta"), py::arg("phi"));

  // Bindings for CoordinateSystem struct
  py::class_<CoordinateSystem>(m, "GridSettings",
                               R"pbdoc(
        Specifying the grid on which the Telescope can request
        a gridded_response or undersampled_response.

        NOTE: the GridSettings class wraps the C++ CoordinateSystem class

       )pbdoc")
      .def(py::init<>(),
           R"pbdoc(
        Initialize CoordinateSystem

        Parameters
        ----------
       )pbdoc")
      .def_readwrite("width", &CoordinateSystem::width,
                     R"pbdoc(
        int: Width of grid in pixels.
       )pbdoc")
      .def_readwrite("height", &CoordinateSystem::height,
                     R"pbdoc(
        int: Height of grid in pixels.
       )pbdoc")
      .def_readwrite("ra", &CoordinateSystem::ra,
                     R"pbdoc(
        double: Right ascension (pointing direction) of telescope [rad].
       )pbdoc")
      .def_readwrite("dec", &CoordinateSystem::dec,
                     R"pbdoc(
        double: Declination (pointing direction) of telescope [rad].
       )pbdoc")
      .def_readwrite("dl", &CoordinateSystem::dl,
                     R"pbdoc(
        double: Grid spacing in RA direction [rad],
        dl is the direction cosine of the delta right ascension
       )pbdoc")
      .def_readwrite("dm", &CoordinateSystem::dm,
                     R"pbdoc(
        double: Grid spacing in Dec direction [rad], where
        dm is the direction cosine of the delta declination.
       )pbdoc")
      .def_readwrite("l_shift", &CoordinateSystem::phase_centre_dl,
                     R"pbdoc(
        double: Shift in l from pointing direction to grid center [rad].
       )pbdoc")
      .def_readwrite("m_shift", &CoordinateSystem::phase_centre_dm,
                     R"pbdoc(
        double: Shift in m from pointing direction to grid center [rad].
       )pbdoc");

  // Bindings for ElementResponseModel enum
  py::enum_<ElementResponseModel>(m, "ElementResponseModel", py::arithmetic(),
                                  "Element Response Model enumeration")
      .value("hamaker", ElementResponseModel::kHamaker,
             R"pbdoc(
        Hamaker element response model
       )pbdoc")
      .value("lobes", ElementResponseModel::kLOBES,
             R"pbdoc(
        LOBEs element response model
       )pbdoc")
      .value("oskar_dipole", ElementResponseModel::kOSKARDipole,
             R"pbdoc(
        SKA dipole element response model
       )pbdoc")
      .value("skala40_spherical", ElementResponseModel::kOSKARSphericalWave,
             R"pbdoc(
        Use SKALA 4.0 element response model. Please note that this response model is somewhat misleadingly named
        OSKAR spherical wave internally (ElementResponseModel::kOSKARSphericalWave). This will be refactored in the future.
       )pbdoc")
      .export_values();

  // Bindings for Options struct
  py::class_<Options>(m, "Options",
                      R"pbdoc(
        Class for specifying some beam forming options.
       )pbdoc")
      .def(py::init<>(),
           R"pbdoc(
        Initialize Options struct

        Parameters
        ----------
       )pbdoc")
      .def_readwrite("coeff_path", &Options::coeff_path,
                     R"pbdoc(
        str: Full path to coefficients file (MWA) or path to directory where
        coefficient files can be found (LOFAR LOBEs model)
       )pbdoc")
      // [TODO] MWA specific
      // .def_readwrite("frequency_interpolation",
      // &Options::frequency_interpolation, "Use frequency interpolation (MWA
      // specific)")
      .def_readwrite("use_differential_beam", &Options::use_differential_beam,
                     R"pbdoc(
        bool: Use differential beam?
       )pbdoc")
      .def_readwrite("use_channel_frequency", &Options::use_channel_frequency,
                     R"pbdoc(
        bool: Use channel frequency?
       )pbdoc")
      .def_readwrite("data_column_name", &Options::data_column_name,
                     R"pbdoc(
        str: Data column name
       )pbdoc")
      .def_readwrite("element_response_model", &Options::element_response_model,
                     R"pbdoc(
        ElementResponseModel: Element response model.
       )pbdoc");
}