#include <pybind11/pybind11.h>

#include "elementresponse.h"
#include "options.h"
#include "coords/coordutils.h"

namespace py = pybind11;
using everybeam::ElementResponseModel;
using everybeam::Options;
using everybeam::coords::CoordinateSystem;

void init_utils(py::module &m) {
  // Bindings for CoordinateSystem struct
  py::class_<CoordinateSystem>(m, "CoordinateSystem",
                               R"pbdoc(
        Class for specifying the dimensions of the grid.
        Will be useful for the gridded response. 
       )pbdoc")
      .def(py::init<>(),
           R"pbdoc(
        Initialize CoordinateSystem
        
        Parameters
        ---------- 
       )pbdoc")
      .def_readwrite("width", &CoordinateSystem::width,
                     R"pbdoc(
        Width of grid in pixels: int
       )pbdoc")
      .def_readwrite("height", &CoordinateSystem::height,
                     R"pbdoc(
        Height of grid in pixels: int
       )pbdoc")
      .def_readwrite("ra", &CoordinateSystem::ra,
                     R"pbdoc(
        Right ascension (pointing direction) of telescope [rad]: double
       )pbdoc")
      .def_readwrite("dec", &CoordinateSystem::dec,
                     R"pbdoc(
        Declination (pointing direction) of telescope [rad]: double
       )pbdoc")
      .def_readwrite("dl", &CoordinateSystem::dl,
                     R"pbdoc(
        Grid size in l direction [rad]: double
       )pbdoc")
      .def_readwrite("dm", &CoordinateSystem::dm,
                     R"pbdoc(
        Grid size in m direction [rad]: double
       )pbdoc")
      .def_readwrite("l_shift", &CoordinateSystem::phase_centre_dl,
                     R"pbdoc(
        Shift in l from pointing direction to grid center [rad]: double
       )pbdoc")
      .def_readwrite("m_shift", &CoordinateSystem::phase_centre_dm,
                     R"pbdoc(
        Shift in m from pointing direction to grid center [rad]: double
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
      .value("osker_spherical", ElementResponseModel::kOSKARSphericalWave,
             R"pbdoc(
        SKA spherical wave element response model
       )pbdoc")
      .export_values();

  // Bindings for Options struct
  py::class_<Options>(m, "Options",
                      R"pbdoc(
        Class for specifying some beam forming options.
       )pbdoc")
      .def(py::init<>())
      // [TODO] MWA specific
      // .def_readwrite("coeff_path", &Options::coeff_path, "Coefficient path to
      // (MWA) coefficient file") .def_readwrite("frequency_interpolation",
      // &Options::frequency_interpolation, "Use frequency interpolation (MWA
      // specific)")
      .def_readwrite("use_differential_beam", &Options::use_differential_beam,
                     R"pbdoc(
        Use differential beam: bool
       )pbdoc")
      .def_readwrite("use_channel_frequency", &Options::use_channel_frequency,
                     R"pbdoc(
        Use channel frequency: bool
       )pbdoc")
      .def_readwrite("data_column_name", &Options::data_column_name,
                     R"pbdoc(
        Data column name: str
       )pbdoc")
      .def_readwrite("element_response_model", &Options::element_response_model,
                     R"pbdoc(
        Element response model: ElementResponseModel
       )pbdoc");
}