// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <array>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "elementresponse.h"
#include "options.h"

namespace py = pybind11;

using everybeam::ElementResponse;
using everybeam::ElementResponseModel;

namespace {

template <typename T, typename U>
void CheckArrayShapes(const py::array_t<T>& left, const py::array_t<U>& right) {
  if (left.ndim() != right.ndim()) {
    throw std::runtime_error("Array dimensions do not match");
  }
  for (py::ssize_t i = 0; i < left.ndim(); ++i) {
    if (left.shape(i) != right.shape(i)) {
      throw std::runtime_error("Array shapes do not match");
    }
  }
}

}  // namespace

void InitElementResponse(py::module& m) {
  // Bindings for ElementResponseModel enum
  py::enum_<ElementResponseModel>(m, "ElementResponseModel", py::arithmetic(),
                                  "Element Response Model enumeration")
      .value("hamaker", ElementResponseModel::kHamaker,
             R"pbdoc(
        Hamaker element response model
       )pbdoc")
      .value("hamaker_lba", ElementResponseModel::kHamakerLba,
             R"pbdoc(
        Hamaker LBA element response model
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

  py::class_<ElementResponse, std::shared_ptr<ElementResponse>>(
      m, "ElementResponse")
      .def_static(
          "create",
          [](const ElementResponseModel model, const std::string& name) {
            return ElementResponse::GetInstance(model, name,
                                                everybeam::Options());
          },
          R"pbdoc(
        Creates an element response object.

        Parameters
        ----------
        model: everybeam.ElementResponseModel
            Model type.
        station: string
            Station name, which is required for various model types.

        Returns
        -------
        An ElementResponse object for the given model.
       )pbdoc",
          py::arg("model"), py::arg("name") = "")
      .def_static(
          "create",
          [](const std::string& model, const std::string& name) {
            return ElementResponse::GetInstance(
                everybeam::ElementResponseModelFromString(model), name,
                everybeam::Options());
          },
          R"pbdoc(
        Creates an element response object.

        Parameters
        ----------
        model: string
            Model type name.
        station: string
            Station name, which is required for various model types.

        Returns
        -------
        An ElementResponse object for the given model.
       )pbdoc",
          py::return_value_policy::reference, py::arg("model"),
          py::arg("name") = "")
      .def_property_readonly(
          "model", [](const ElementResponse& self) { return self.GetModel(); },
          R"pbdoc(
        Returns
        -------
        everybeam.ElementResponseModel
            The model of the element response
       )pbdoc")
      .def(
          "response",
          [](const ElementResponse& self, const py::array_t<double>& freq,
             const py::array_t<double>& theta, const py::array_t<double>& phi) {
            CheckArrayShapes(freq, theta);
            CheckArrayShapes(freq, phi);
            std::vector<py::ssize_t> result_shape;
            for (py::ssize_t i = 0; i < freq.ndim(); ++i) {
              result_shape.push_back(freq.shape(i));
            }
            result_shape.push_back(2);
            result_shape.push_back(2);
            py::array_t<std::complex<double>> result(result_shape);

            const double* freq_data = freq.data();
            const double* theta_data = theta.data();
            const double* phi_data = phi.data();
            auto result_data =
                reinterpret_cast<aocommon::MC2x2*>(result.mutable_data());
            // Verify that the reinterpret_cast is allowed.
            static_assert(sizeof(aocommon::MC2x2) == 2 * 4 * sizeof(double));

            for (py::ssize_t i = 0; i < freq.size(); ++i) {
              result_data[i] =
                  self.Response(freq_data[i], theta_data[i], phi_data[i]);
            }
            return result;
          },
          R"pbdoc(
        Compute the response.

        This function supports array arguments. When using array arguments,
        all arguments must have the same array shape. The return value is
        an array with two extra dimensions, and contains the Jones matrix
        for each input array element.

        Parameters
        ----------
        freq: double
            Frequency of the plane wave (Hz).
        theta: double
            Angle wrt. z-axis (rad)
        phi: double
            Angle in the xy-plane wrt. x-axis  (rad)

        Returns
        -------
        np.ndarray, np.complex64
            2-D numpy array ``[2, 2]`` containing the Jones matrix.
       )pbdoc",
          py::arg("freq"), py::arg("theta"), py::arg("phi"))
      .def(
          "response",
          [](const ElementResponse& self, const py::array_t<int>& element_id,
             const py::array_t<double>& freq, const py::array_t<double>& theta,
             const py::array_t<double>& phi) {
            CheckArrayShapes(element_id, freq);
            CheckArrayShapes(freq, theta);
            CheckArrayShapes(freq, phi);
            std::vector<py::ssize_t> result_shape;
            for (py::ssize_t i = 0; i < freq.ndim(); ++i) {
              result_shape.push_back(freq.shape(i));
            }
            result_shape.push_back(2);
            result_shape.push_back(2);
            py::array_t<std::complex<double>> result(result_shape);

            const int* id_data = element_id.data();
            const double* freq_data = freq.data();
            const double* theta_data = theta.data();
            const double* phi_data = phi.data();
            auto result_data =
                reinterpret_cast<aocommon::MC2x2*>(result.mutable_data());

            for (py::ssize_t i = 0; i < freq.size(); ++i) {
              result_data[i] = self.Response(id_data[i], freq_data[i],
                                             theta_data[i], phi_data[i]);
            }
            return result;
          },
          R"pbdoc(
        Compute the response for a specific element.

        This function supports array arguments. When using array arguments,
        all arguments must have the same array shape. The return value is
        an array with two extra dimensions, and contains the Jones matrix
        for each input array element.

        Parameters
        ----------
        element_id: int
            ID of an element.
        freq: double
            Frequency of the plane wave (Hz).
        theta: double
            Angle wrt. z-axis (rad)
        phi: double
            Angle in the xy-plane wrt. x-axis  (rad)

        Returns
        -------
        np.ndarray, np.complex64
            2-D numpy array ``[2, 2]`` containing the Jones matrix.
       )pbdoc",
          py::arg("element_id"), py::arg("freq"), py::arg("theta"),
          py::arg("phi"));
}
