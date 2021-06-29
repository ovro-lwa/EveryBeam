// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "load.h"
#include "station.h"
#include "griddedresponse/griddedresponse.h"
#include "common/types.h"
#include "coords/coordutils.h"
#include "coords/itrfconverter.h"
#include "coords/itrfdirection.h"
#include "telescope/lofar.h"
#include "telescope/phasedarray.h"

#include <aocommon/matrix2x2.h>

namespace py = pybind11;

using casacore::MeasurementSet;

using everybeam::Options;
using everybeam::Station;
using everybeam::vector3r_t;
using everybeam::coords::CoordinateSystem;
using everybeam::coords::ITRFConverter;
using everybeam::coords::SetITRFVector;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::telescope::LOFAR;
using everybeam::telescope::PhasedArray;
using everybeam::telescope::Telescope;

namespace {
// Convert pyarray of size 3 to vector3r_t
vector3r_t np2vector3r_t(const py::array_t<double> &pyarray) {
  auto r = pyarray.unchecked<1>();
  if (r.size() != 3) {
    throw std::runtime_error("Pyarry is of incorrect size, must be 3.");
  }
  return vector3r_t{r[0], r[1], r[2]};
}

// Cast aocommon::MC2x2 to py::array
py::array_t<std::complex<double>> cast_matrix(const aocommon::MC2x2 &matrix) {
  // Solution from: https://github.com/pybind/pybind11/issues/1299
  // Reinterpret cast is needed to "flatten" the nested std::array
  auto mat_ptr = reinterpret_cast<const std::complex<double> *>(matrix.Data());
  return py::array_t<std::complex<double>>(std::vector<ptrdiff_t>{2, 2},
                                           mat_ptr);
}

// Cast vector of aocommon::MC2x2 to numpy tensor
py::array_t<std::complex<double>> cast_tensor(
    const std::vector<aocommon::MC2x2> &matrix,
    const std::vector<size_t> &layout) {
  size_t total_size = 1;
  for (size_t rank_size : layout) {
    total_size *= rank_size;
  }

  if (total_size != (matrix.size() * 4)) {
    throw std::runtime_error("Casting mismatching shapes");
  }
  // Reinterpret_cast needed to flatten the nested std::array
  auto mat_ptr = reinterpret_cast<const std::complex<double> *>(matrix.data());
  std::vector<ptrdiff_t> np_layout(layout.begin(), layout.end());
  return py::array_t<std::complex<double>>(np_layout, mat_ptr);
}

// Convenience method throwing an error if provided station
// index exceeds number of stations
void check_station_index(size_t idx, size_t idx_max,
                         const std::string &prefix) {
  if (idx >= idx_max) {
    throw std::runtime_error(
        prefix + ": Requested station index exceeds number of stations.");
  }
}

// Convenience function for applying the differential beam response (lhs) to the
// response (rhs)
void apply_differential_beam(aocommon::MC2x2 &lhs, const aocommon::MC2x2 &rhs) {
  // Computes lhs = lhs^{-1} * rhs by
  // 1) computing the inverse of response (in preapplied beam direction), if it
  // exists...
  if (!lhs.Invert()) {
    lhs = aocommon::MC2x2::Zero();
  }
  // 2) and multiplying the inverse with response (rhs)
  lhs *= rhs;
}

}  // namespace

/**
 * @brief Trampoline class for Telescope. Needed for correct overloadin
 * in python. Not yet needed, but can become useful in the near future.
 * See
 * https://pybind11.readthedocs.io/en/stable/advanced/classes.html?highlight=trampoline#classes
 *
 */
class PyTelescope : public Telescope {
  /* Inherit the constructors */
  using Telescope::Telescope;

  GriddedResponse *GetGriddedResponseWrapper(
      const CoordinateSystem &coordinate_system) {
    PYBIND11_OVERLOAD_PURE(GriddedResponse *, Telescope, GetGriddedResponse,
                           coordinate_system);
  }
};

std::unique_ptr<LOFAR> create_lofar(const std::string &name,
                                    const everybeam::Options &options) {
  MeasurementSet ms(name);

  if (everybeam::GetTelescopeType(ms) !=
      everybeam::TelescopeType::kLofarTelescope) {
    throw std::runtime_error(
        "LOFAR telescope requested, but specified path name does not contain a "
        "LOFAR MS.");
  }
  std::unique_ptr<LOFAR> telescope =
      std::unique_ptr<LOFAR>(new LOFAR(ms, options));
  return telescope;
}

void init_telescope(py::module &m) {
  py::class_<Telescope, PyTelescope>(m, "Telescope")
      .def_property_readonly("is_time_relevant", &Telescope::GetIsTimeRelevant,
                             R"pbdoc(
        Is time relevant for this telescope?

        Returns
        -------
        bool
       )pbdoc")
      .def_property_readonly("nr_stations", &Telescope::GetNrStations,
                             R"pbdoc(
        Retrieve the number of stations.

        Returns
        -------
        int
       )pbdoc")
      .def_property_readonly("options", &Telescope::GetOptions,
                             R"pbdoc(
        Retrieve the specified options.

        Returns
        -------
        pyeverybeam.Options
       )pbdoc");

  py::class_<PhasedArray, Telescope>(m, "PhasedArray")
      .def_property_readonly("nr_channels", &PhasedArray::GetNrChannels,
                             R"pbdoc(
        Retrieve the number of channels.

        Returns
        -------
        int
       )pbdoc")
      .def("station_name",
           [](PhasedArray &self, size_t idx) {
             if (idx >= self.GetNrStations()) {
               throw std::runtime_error(
                   "Requested station index exceeds number of stations.");
             }
             const Station &station =
                 static_cast<const Station &>(*(self.GetStation(idx).get()));
             return station.GetName();
           })
      .def("channel_frequency", &PhasedArray::GetChannelFrequency,
           R"pbdoc(
        Retrieve channel frequency for a given (zero-based) channel index.

        Parameters
        ----------
        channel_index: int
            Channel index

        Returns
        -------
        float
       )pbdoc",
           py::arg("channel_index"))
      // Corresponds to evaluate0 in lofarbeam
      .def(
          "station_response",
          [](PhasedArray &self, double time,
             bool rotate) -> py::array_t<std::complex<double>> {
            vector3r_t direction, station0, tile0;
            ITRFConverter itrf_converter(time);
            SetITRFVector(
                itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                direction);
            SetITRFVector(itrf_converter.ToDirection(self.GetDelayDirection()),
                          station0);
            SetITRFVector(
                itrf_converter.ToDirection(self.GetTileBeamDirection()), tile0);

            const size_t nr_stations = self.GetNrStations();
            const size_t nr_channels = self.GetNrChannels();
            std::vector<aocommon::MC2x2> response(nr_stations * nr_channels);
            for (size_t station_idx = 0; station_idx < nr_stations;
                 ++station_idx) {
              const Station &station = static_cast<const Station &>(
                  *(self.GetStation(station_idx).get()));
              for (size_t channel_idx = 0; channel_idx < nr_channels;
                   ++channel_idx) {
                const size_t idx = station_idx * nr_channels + channel_idx;
                const double freq = self.GetChannelFrequency(channel_idx);
                const double freq0 = self.GetOptions().use_channel_frequency
                                         ? freq
                                         : self.GetSubbandFrequency();

                if (self.GetOptions().use_differential_beam) {
                  // Exploiting: R^{-1}R = I
                  response[idx] = aocommon::MC2x2::Unity();
                } else {
                  response[idx] = station.Response(time, freq, direction, freq0,
                                                   station0, tile0, rotate);
                }
              }
            }
            return py::array_t<std::complex<double>>{cast_tensor(
                response, std::vector<size_t>{nr_stations, nr_channels, 2, 2})};
          },
          R"pbdoc(
        Get station response in beam former direction for all stations and channels.

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        rotate: bool, optional
            Apply paralactic angle rotation? [True/False] Defaults to True

        Returns
        -------
        np.ndarray
            rank 4 numpy array of shape (nr_stations, nr_channels, 2, 2)
       )pbdoc",
          py::arg("time"), py::arg("rotate") = true)
      // Corresponds to evaluate1 in lofarbeam
      .def(
          "station_response",
          [](PhasedArray &self, double time, size_t idx,
             bool rotate) -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "station_response");
            vector3r_t direction, station0, tile0;
            ITRFConverter itrf_converter(time);
            SetITRFVector(
                itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                direction);
            SetITRFVector(itrf_converter.ToDirection(self.GetDelayDirection()),
                          station0);
            SetITRFVector(
                itrf_converter.ToDirection(self.GetTileBeamDirection()), tile0);

            const Station &station =
                static_cast<const Station &>(*(self.GetStation(idx).get()));

            const size_t nr_channels = self.GetNrChannels();
            std::vector<aocommon::MC2x2> response(nr_channels);
            for (size_t idx = 0; idx < nr_channels; ++idx) {
              const double freq = self.GetChannelFrequency(idx);
              const double freq0 = self.GetOptions().use_channel_frequency
                                       ? freq
                                       : self.GetSubbandFrequency();
              if (self.GetOptions().use_differential_beam) {
                // Exploiting R^{-1}R = I
                response[idx] = aocommon::MC2x2::Unity();
              } else {
                response[idx] = station.Response(time, freq, direction, freq0,
                                                 station0, tile0, rotate);
              }
            }
            return py::array_t<std::complex<double>>{
                cast_tensor(response, std::vector<size_t>{nr_channels, 2, 2})};
          },
          R"pbdoc(
        Get station response in beam former direction for all channels.

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        station_idx: int
            Get response for station index
        rotate: bool, optional
            Apply paralactic angle rotation? [True/False] Defaults to True

        Returns
        -------
        np.ndarray
            rank 3 numpy array of shape (nr_channels, 2, 2)
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("rotate") = true)
      // Corresponds to evaluate2 in lofarbeam
      .def(
          "station_response",
          [](PhasedArray &self, double time, size_t idx, size_t channel,
             bool rotate) -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "station_response");
            if (channel >= self.GetNrChannels()) {
              throw std::runtime_error(
                  "Requested channel index exceeds channel count.");
            }

            const double freq = self.GetChannelFrequency(channel);
            const double freq0 = self.GetOptions().use_channel_frequency
                                     ? freq
                                     : self.GetSubbandFrequency();

            vector3r_t direction, station0, tile0;
            ITRFConverter itrf_converter(time);
            SetITRFVector(
                itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                direction);
            SetITRFVector(itrf_converter.ToDirection(self.GetDelayDirection()),
                          station0);
            SetITRFVector(
                itrf_converter.ToDirection(self.GetTileBeamDirection()), tile0);

            const Station &station =
                static_cast<const Station &>(*(self.GetStation(idx).get()));

            const aocommon::MC2x2 response =
                self.GetOptions().use_differential_beam
                    ? aocommon::MC2x2::Unity()
                    : station.Response(time, freq, direction, freq0, station0,
                                       tile0, rotate);
            return py::array_t<std::complex<double>>{cast_matrix(response)};
          },
          R"pbdoc(
        Get station response in beam former direction for specified channel.

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        station_idx: int
            Get response for station index
        channel_idx: int
            Index of channel.
        rotate: bool, optional
            Apply paralactic angle rotation? [True/False] Defaults to True

        Returns
        -------
        np.ndarray
            Response (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("channel_idx"),
          py::arg("rotate") = true)
      // Corresponds to evaluate3 in lofarbeam
      .def(
          "station_response",
          [](PhasedArray &self, double time, size_t idx, double freq,
             bool rotate) -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "station_response");
            const double freq0 = self.GetOptions().use_channel_frequency
                                     ? freq
                                     : self.GetSubbandFrequency();

            vector3r_t direction, station0, tile0;
            ITRFConverter itrf_converter(time);
            SetITRFVector(
                itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                direction);
            SetITRFVector(itrf_converter.ToDirection(self.GetDelayDirection()),
                          station0);
            SetITRFVector(
                itrf_converter.ToDirection(self.GetTileBeamDirection()), tile0);

            const Station &station =
                static_cast<const Station &>(*(self.GetStation(idx).get()));

            const aocommon::MC2x2 response =
                self.GetOptions().use_differential_beam
                    ? aocommon::MC2x2::Unity()
                    : station.Response(time, freq, direction, freq0, station0,
                                       tile0, rotate);
            return py::array_t<std::complex<double>>{cast_matrix(response)};
          },
          R"pbdoc(
        Get station response in beam former direction for specified frequency.

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        station_idx: int
            Get response for station index
        freq: float
            Frequency of the plane wave (Hz)
        rotate: bool, optional
            Apply paralactic angle rotation? [True/False] Defaults to True

        Returns
        -------
        np.ndarray
            Response (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("freq"),
          py::arg("rotate") = true)
      // Corresponds to evaluate4 from lofarbeam
      .def(
          "station_response",
          [](PhasedArray &self, double time, size_t idx, double freq,
             const py::array_t<double> pydirection,
             const py::array_t<double> pystation0,
             const py::array_t<double> pytile0,
             bool rotate) -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "station_response");
            const vector3r_t direction = np2vector3r_t(pydirection);
            const vector3r_t station0 = np2vector3r_t(pystation0);
            const vector3r_t tile0 = np2vector3r_t(pytile0);
            const double freq0 = self.GetOptions().use_channel_frequency
                                     ? freq
                                     : self.GetSubbandFrequency();

            const Station &station =
                static_cast<const Station &>(*(self.GetStation(idx).get()));
            const aocommon::MC2x2 response = station.Response(
                time, freq, direction, freq0, station0, tile0, rotate);

            if (self.GetOptions().use_differential_beam) {
              vector3r_t diff_beam_centre;
              ITRFConverter itrf_converter(time);
              SetITRFVector(
                  itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                  diff_beam_centre);

              aocommon::MC2x2 response_diff_beam = station.Response(
                  time, freq, diff_beam_centre, freq0, station0, tile0, rotate);
              apply_differential_beam(response_diff_beam, response);

              return py::array_t<std::complex<double>>{
                  cast_matrix(response_diff_beam)};
            } else {
              return py::array_t<std::complex<double>>{cast_matrix(response)};
            }
          },
          R"pbdoc(
        Get station response in user-specified direction

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        station_idx: int
            Get response for station index
        freq: float
            Frequency of the plane wave (Hz)
        direction: np.1darray
            Direction of arrival (ITRF, m)
        station0: np.1darray
            Station beam former reference direction (ITRF, m)
        tile0: np.1darray
            Tile beam former reference direction (ITRF, m)
        rotate: bool, optional
            Apply paralactic angle rotation? [True/False] Defaults to True

        Returns
        -------
        np.ndarray
            Response (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("freq"),
          py::arg("direction"), py::arg("station0_direction"),
          py::arg("tile0_direction"), py::arg("rotate") = true)
      // Same as previous, but station0 direction and
      // tile0 direction coincide.
      .def(
          "station_response",
          [](PhasedArray &self, double time, size_t idx, double freq,
             const py::array_t<double> pydirection,
             const py::array_t<double> pystation0,
             bool rotate) -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "station_response");
            const vector3r_t direction = np2vector3r_t(pydirection);
            const vector3r_t station0 = np2vector3r_t(pystation0);
            // tile0 direction identical to station0

            const double freq0 = self.GetOptions().use_channel_frequency
                                     ? freq
                                     : self.GetSubbandFrequency();

            const Station &station =
                static_cast<const Station &>(*(self.GetStation(idx).get()));
            const aocommon::MC2x2 response = station.Response(
                time, freq, direction, freq0, station0, station0, rotate);

            if (self.GetOptions().use_differential_beam) {
              vector3r_t diff_beam_centre;
              ITRFConverter itrf_converter(time);
              SetITRFVector(
                  itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                  diff_beam_centre);

              aocommon::MC2x2 response_diff_beam =
                  station.Response(time, freq, diff_beam_centre, freq0,
                                   station0, station0, rotate);

              apply_differential_beam(response_diff_beam, response);
              return py::array_t<std::complex<double>>{
                  cast_matrix(response_diff_beam)};
            } else {
              return py::array_t<std::complex<double>>{cast_matrix(response)};
            }
          },
          R"pbdoc(
        Get station response in user-specified direction

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        station_idx: int
            station index
        freq: float
            Frequency of the plane wave (Hz)
        direction: np.1darray
            Direction of arrival (ITRF, m)
        station0: np.1darray
            Station beam former reference direction (ITRF, m)
        rotate: bool, optional
            Apply paralactic angle rotation? [True/False] Defaults to True

        Returns
        -------
        np.ndarray
            Response (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("freq"),
          py::arg("direction"), py::arg("station0_direction"),
          py::arg("rotate") = true)
      .def(
          "element_response",
          [](PhasedArray &self, double time, size_t idx, size_t element_idx,
             double freq, const py::array_t<double> pydirection, bool is_local,
             bool rotate) -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "element_response");
            // TODO: need such a check for element index too

            const vector3r_t direction = np2vector3r_t(pydirection);
            const Station &station =
                static_cast<const Station &>(*(self.GetStation(idx).get()));
            const aocommon::MC2x2 response = station.ComputeElementResponse(
                time, freq, direction, element_idx, is_local, rotate);

            if (self.GetOptions().use_differential_beam) {
              vector3r_t diff_beam_centre;
              ITRFConverter itrf_converter(time);
              SetITRFVector(
                  itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                  diff_beam_centre);

              aocommon::MC2x2 response_diff_beam =
                  station.ComputeElementResponse(time, freq, diff_beam_centre,
                                                 element_idx, is_local, rotate);
              apply_differential_beam(response_diff_beam, response);
              return py::array_t<std::complex<double>>{
                  cast_matrix(response_diff_beam)};
            } else {
              return py::array_t<std::complex<double>>{cast_matrix(response)};
            }
          },
          R"pbdoc(
        Get element response given a station and an element in prescribed direction.

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        station_idx: int
            station index
        element_idx: int
            element index
        freq: float
            Frequency of the plane wave (Hz)
        direction: np.1darray
            Direction of arrival either in ITRF (m) or local East-North-Up (m)
        is_local: bool, optional
            Is the specified direction in local East-North-Up? If not, global coordinate
            system is assumed. [True/False] Defaults to False.
        rotate: bool, optional
            Apply paralactic angle rotation? [True/False] Defaults to True

        Returns
        -------
        np.ndarray
            Response (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("element_idx"),
          py::arg("freq"), py::arg("direction"), py::arg("is_local") = false,
          py::arg("rotate") = true)
      .def(
          "element_response",
          [](PhasedArray &self, double time, size_t idx, double freq,
             const py::array_t<double> pydirection, bool is_local,
             bool rotate) -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "element_response");
            // TODO: need such a check for element index too

            const vector3r_t direction = np2vector3r_t(pydirection);
            const Station &station =
                static_cast<const Station &>(*(self.GetStation(idx).get()));
            const aocommon::MC2x2 response = station.ComputeElementResponse(
                time, freq, direction, is_local, rotate);

            if (self.GetOptions().use_differential_beam) {
              vector3r_t diff_beam_centre;
              ITRFConverter itrf_converter(time);
              SetITRFVector(
                  itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                  diff_beam_centre);

              aocommon::MC2x2 response_diff_beam =
                  station.ComputeElementResponse(time, freq, diff_beam_centre,
                                                 is_local, rotate);
              apply_differential_beam(response_diff_beam, response);
              return py::array_t<std::complex<double>>{
                  cast_matrix(response_diff_beam)};
            } else {
              return py::array_t<std::complex<double>>{cast_matrix(response)};
            }

            // return py::array_t<std::complex<double>>{cast_matrix(response)};
          },
          R"pbdoc(
        Get element response given a station index, in prescribed direction.

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        station_idx: int
            station index
        freq: float
            Frequency of the plane wave (Hz)
        direction: np.1darray
            Direction of arrival either in ITRF (m) or local East-North-Up (m)
        is_local: bool, optional
            Is the specified direction in local East-North-Up? If not, global coordinate
            system is assumed. [True/False] Defaults to False.
        rotate: bool, optional
            Apply paralactic angle rotation? [True/False] Defaults to True

        Returns
        -------
        np.ndarray
            Response (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("freq"),
          py::arg("direction"), py::arg("is_local") = false,
          py::arg("rotate") = true)
      .def(
          "array_factor",
          [](PhasedArray &self, double time, size_t idx, double freq,
             const py::array_t<double> pydirection,
             const py::array_t<double> pystation0)
              -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "array_factor");
            const vector3r_t direction = np2vector3r_t(pydirection);
            const vector3r_t station0 = np2vector3r_t(pystation0);
            const double freq0 = self.GetOptions().use_channel_frequency
                                     ? freq
                                     : self.GetSubbandFrequency();

            const Station &station =
                static_cast<const Station &>(*(self.GetStation(idx).get()));
            const aocommon::MC2x2Diag response_diag = station.ArrayFactor(
                time, freq, direction, freq0, station0, station0);
            // Diagonal to 2x2 matrix
            const aocommon::MC2x2 response(response_diag[0], 0.0, 0.0,
                                           response_diag[1]);
            return py::array_t<std::complex<double>>{cast_matrix(response)};
          },
          R"pbdoc(
        Get array factor for a given station in prescribed direction.

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        station_idx: int
            station index
        freq: float
            Frequency of the plane wave (Hz)
        direction: np.1darray
            Direction of arrival in ITRF (m)
        station0_direction: np.1darray
            Station beam former reference direction (ITRF, m)

        Returns
        -------
        np.ndarray
            Response diagonal (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("freq"),
          py::arg("direction"), py::arg("station0_direction"))
      .def(
          "array_factor",
          [](PhasedArray &self, double time, size_t idx, double freq,
             const py::array_t<double> pydirection,
             const py::array_t<double> pystation0,
             const py::array_t<double> pytile0)
              -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "array_factor");
            const vector3r_t direction = np2vector3r_t(pydirection);
            const vector3r_t station0 = np2vector3r_t(pystation0);
            const vector3r_t tile0 = np2vector3r_t(pytile0);
            const double freq0 = self.GetOptions().use_channel_frequency
                                     ? freq
                                     : self.GetSubbandFrequency();

            const Station &station =
                static_cast<const Station &>(*(self.GetStation(idx).get()));
            const aocommon::MC2x2Diag response_diag = station.ArrayFactor(
                time, freq, direction, freq0, station0, tile0);
            // Diagonal to 2x2 matrix
            const aocommon::MC2x2 response(response_diag[0], 0.0, 0.0,
                                           response_diag[1]);
            return py::array_t<std::complex<double>>{cast_matrix(response)};
          },
          R"pbdoc(
        Get array factor for a given station in prescribed direction.

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        station_idx: int
            station index
        freq: float
            Frequency of the plane wave (Hz)
        direction: np.1darray
            Direction of arrival in ITRF (m)
        station0_direction: np.1darray
            Station beam former reference direction (ITRF, m)
        tile0_direction: np.1darray
            Tile beam former reference direction (ITRF, m)

        Returns
        -------
        np.ndarray
            Response diagonal (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("freq"),
          py::arg("direction"), py::arg("station0_direction"),
          py::arg("tile0_direction"));

  py::class_<LOFAR, PhasedArray>(m, "LOFAR").def(py::init(&create_lofar));

  // TODO: other telescopes:
  // py::class_<OSKAR, PhasedArray>(m, "OSKAR");

  // py::class_<MWA, Telescope>(m, "MWA");

  // py::class_<Dish, Telescope>(m, "Dish");

  // py::class_<VLA, Dish>(m, "MWA");
}
