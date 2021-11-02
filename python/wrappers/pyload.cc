// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>
#include <memory>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "load.h"
#include "options.h"

using casacore::MeasurementSet;

using everybeam::GetTelescopeType;
using everybeam::Load;
using everybeam::telescope::Telescope;

namespace py = pybind11;

// Wrapper around the everybeam::Load method
std::unique_ptr<Telescope> pyload_telescope(
    const std::string& name, const std::string& data_column = "DATA",
    bool differential_beam = false, bool channel_frequency = true,
    const std::string& element_response_model = "hamaker",
    const std::string& coeff_path = "") {
  // Load measurement set
  MeasurementSet ms(name);

  if (everybeam::GetTelescopeType(ms) !=
          everybeam::TelescopeType::kLofarTelescope &&
      everybeam::GetTelescopeType(ms) !=
          everybeam::TelescopeType::kOSKARTelescope) {
    throw std::runtime_error(
        "Currently the python bindings only support the LOFAR and the OSKAR "
        "telescope");
  }

  // Fill everybeam options
  everybeam::Options options;

  // LOFAR related
  std::string element_response_tmp = element_response_model;
  std::for_each(element_response_tmp.begin(), element_response_tmp.end(),
                [](char& c) { c = ::toupper(c); });
  everybeam::ElementResponseModel element_response_enum;
  if (element_response_tmp == "HAMAKER")
    element_response_enum = everybeam::ElementResponseModel::kHamaker;
  else if (element_response_tmp == "LOBES")
    element_response_enum = everybeam::ElementResponseModel::kLOBES;
  else if (element_response_tmp == "OSKAR_DIPOLE")
    element_response_enum = everybeam::ElementResponseModel::kOSKARDipole;
  else if (element_response_tmp == "SKALA40_WAVE")
    element_response_enum =
        everybeam::ElementResponseModel::kOSKARSphericalWave;
  else {
    std::stringstream message;
    message << "The specified element response model " << element_response_model
            << " is not implemented.";
    throw std::runtime_error(message.str());
  }
  options.data_column_name = data_column;
  options.use_differential_beam = differential_beam;
  options.use_channel_frequency = channel_frequency;
  options.element_response_model = element_response_enum;
  options.coeff_path = coeff_path;

  // Return the telescope
  std::unique_ptr<Telescope> telescope = Load(ms, options);
  return telescope;
}

void init_load(py::module& m) {
  m.def("load_telescope", &pyload_telescope, R"pbdoc(
        Load telescope from measurement set (MS)

        Parameters
        ----------
        name: str
            Path to MS
        data_column: str, optional
            Data column that should
        use_differential_beam: bool, optional
            Use differential beam? Defaults to False
        use_channel_frequency: bool, optional
            Use channel frequency? Defaults to True.
        element_response_model: str
            Specify the element response model, should be any of
            ["hamaker", "lobes", "oskar_dipole", "skala40_wave"]
            Please note that the SKALA40 Wave model is
            currently named OSKAR Spherical Wave in the EveryBeam internals.
            This will be refactored to SKALA40_WAVE in the future.

        Returns
        -------
        Telescope object
       )pbdoc",
        py::arg("name"), py::arg("data_column") = "DATA",
        py::arg("use_differential_beam") = false,
        py::arg("use_channel_frequency") = true,
        py::arg("element_response_model") = "hamaker",
        py::arg("coeff_path") = "");
}