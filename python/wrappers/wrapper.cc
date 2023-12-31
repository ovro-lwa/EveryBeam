// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>

namespace py = pybind11;

void InitElementResponse(py::module&);
void init_load(py::module&);
void init_lobes(py::module&);
void init_telescope(py::module&);
void init_utils(py::module&);

PYBIND11_MODULE(everybeam, m) {
  m.doc() = R"pbdoc(
   pyeverybeam provides python-wrappers for the everybeam beam response library
  )pbdoc";
  // InitElementResponse must come before init_utils, since init_utils uses
  // ElementResponseModel for Options.
  InitElementResponse(m);
  // init_utils must come before init_load, since init_load uses
  // BeamNormalisationMode.
  init_utils(m);
  init_load(m);
  init_telescope(m);
  init_lobes(m);
}