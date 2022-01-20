// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>
#include <memory>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "load.h"
#include "options.h"

using casacore::MeasurementSet;

using everybeam::BeamNormalisationMode;
using everybeam::GetTelescopeType;
using everybeam::Load;
using everybeam::telescope::Telescope;

namespace py = pybind11;

// Wrapper around the everybeam::Load method
std::unique_ptr<Telescope> pyload_telescope(
    const std::string& name, const std::string& data_column,
    BeamNormalisationMode beam_normalisation_mode, bool use_channel_frequency,
    const std::string& element_response_model, const std::string& coeff_path) {
  // Load measurement set
  MeasurementSet ms(name);

  const everybeam::TelescopeType telescope_type =
      everybeam::GetTelescopeType(ms);

  switch (telescope_type) {
    case everybeam::TelescopeType::kAARTFAAC:
    case everybeam::TelescopeType::kLofarTelescope:
    case everybeam::TelescopeType::kOSKARTelescope:
    case everybeam::TelescopeType::kSkaMidTelescope:
      break;
    default:
      throw std::runtime_error(
          "Currently the python bindings only support AARTFAAC (LBA), LOFAR, "
          "OSKAR (SKALA40) and SKA-MID observations");
  }

  // Fill everybeam options
  everybeam::Options options;

  std::string element_response_tmp = element_response_model;
  std::for_each(element_response_tmp.begin(), element_response_tmp.end(),
                [](char& c) { c = ::toupper(c); });
  if (element_response_tmp == "HAMAKER") {
    options.element_response_model = everybeam::ElementResponseModel::kHamaker;
  } else if (element_response_tmp == "HAMAKER_LBA") {
    options.element_response_model =
        everybeam::ElementResponseModel::kHamakerLba;
  } else if (element_response_tmp == "LOBES") {
    options.element_response_model = everybeam::ElementResponseModel::kLOBES;
  } else if (element_response_tmp == "OSKAR_DIPOLE") {
    options.element_response_model =
        everybeam::ElementResponseModel::kOSKARDipole;
  } else if (element_response_tmp == "SKALA40_WAVE") {
    options.element_response_model =
        everybeam::ElementResponseModel::kOSKARSphericalWave;
  } else if (element_response_tmp == "SKAMID_ANALYTICAL") {
    options.element_response_model =
        everybeam::ElementResponseModel::kSkaMidAnalytical;
  } else {
    std::stringstream message;
    message << "The specified element response model " << element_response_model
            << " is not implemented.";
    throw std::runtime_error(message.str());
  }
  options.data_column_name = data_column;
  options.beam_normalisation_mode = beam_normalisation_mode;
  options.use_channel_frequency = use_channel_frequency;
  options.coeff_path = coeff_path;

  return Load(ms, options);
}

void init_load(py::module& m) {
  m.def(
      "load_telescope",
      [](const std::string& name, const std::string& data_column,
         bool use_differential_beam, bool use_channel_frequency,
         const std::string& element_response_model,
         const std::string& coeff_path = "") -> std::unique_ptr<Telescope> {
        BeamNormalisationMode beam_normalisation_mode =
            use_differential_beam ? BeamNormalisationMode::kPreApplied
                                  : BeamNormalisationMode::kNone;
        return pyload_telescope(name, data_column, beam_normalisation_mode,
                                use_channel_frequency, element_response_model,
                                coeff_path);
      },
      R"pbdoc(
        Load telescope from measurement set (MS)

        This version has a simple on/off toggle for beam normalisation through
        the use_differential_beam parameter

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

  m.def("load_telescope", &pyload_telescope, R"pbdoc(
        Load telescope from measurement set (MS)

        This version allows more fine grained control over the normalisation
        of the beam through the beam_normalisation_mode parameter.
        (needed by the DP3 python step implemented in idgcaldpstep.py in the IDG library)

        Parameters
        ----------
        name: str
            Path to MS
        data_column: str, optional
            Data column that should
        beam_normalisation_mode : BeamNormalisationMode, optional
            Defaults to BeamNormalisationMode.none (no normalisation)
            see BeamNormalisationMode enum
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
        py::arg("beam_normalisation_mode") = BeamNormalisationMode::kNone,
        py::arg("use_channel_frequency") = true,
        py::arg("element_response_model") = "hamaker",
        py::arg("coeff_path") = "");
}