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

// Cast complex buffer to numpy array
template <typename T>
py::array_t<std::complex<T>> cast_tensor(const std::complex<T> *buffer,
                                         const std::vector<size_t> &layout) {
  std::vector<ptrdiff_t> np_layout(layout.begin(), layout.end());
  return py::array_t<std::complex<T>>(np_layout, buffer);
}

// Cast aocommon::MC2x2 to py::array
py::array_t<std::complex<double>> cast_matrix(const aocommon::MC2x2 &matrix) {
  // Solution from: https://github.com/pybind/pybind11/issues/1299
  return cast_tensor<double>(matrix.Data(), {2, 2});
}

// Cast vector of aocommon complex matrices (MC2x2/MC4x4) to numpy tensor
// with some additional size checks
template <size_t nelem, typename T>
py::array_t<std::complex<double>> cast_tensor_mc(
    const std::vector<T> &matrix, const std::vector<size_t> &layout) {
  size_t total_size = 1;
  for (size_t rank_size : layout) {
    total_size *= rank_size;
  }

  if (total_size != (matrix.size() * nelem)) {
    throw std::runtime_error("Casting mismatching shapes");
  }

  // Reinterpret_cast needed to flatten the nested std::array
  auto buffer = reinterpret_cast<const std::complex<double> *>(matrix.data());
  return cast_tensor<double>(buffer, layout);
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

// Convenience method to convert a buffer of [npixels * 16]
// (where 16 the number of entries of an aocmmon::HMC4x4 matrix
// to a vector of MC4x4 matrices
std::vector<aocommon::MC4x4> hmc_to_mc(double *buffer, size_t npixels) {
  std::vector<aocommon::MC4x4> vec_mc4x4(npixels);
  for (size_t pixel = 0; pixel != npixels; ++pixel) {
    std::array<double, 16> tmp;
    // Collect the components of the Hermitian matrix
    for (size_t i = 0; i != 16; ++i) tmp[i] = buffer[pixel + i * npixels];
    aocommon::HMC4x4 hmc4x4(tmp.data());
    vec_mc4x4[pixel] = hmc4x4.ToMatrix();
  }
  return vec_mc4x4;
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
       )pbdoc")
      .def(
          "gridded_response",
          [](Telescope &self, const CoordinateSystem &coordinate_system,
             double time, double freq, size_t station_idx,
             size_t field_idx) -> py::array_t<std::complex<float>> {
            check_station_index(station_idx, self.GetNrStations(),
                                "gridded_response");
            std::unique_ptr<GriddedResponse> grid_response =
                self.GetGriddedResponse(coordinate_system);
            std::vector<std::complex<float>> buffer(
                grid_response->GetStationBufferSize(1));
            grid_response->FullResponse(buffer.data(), time, freq, station_idx,
                                        field_idx);
            return cast_tensor<float>(
                buffer.data(),
                {coordinate_system.height, coordinate_system.width, 2, 2});
          },
          R"pbdoc(
        Compute the gridded response for a single station.

        Parameters
        ----------
        coordinate_system: everybeam.CoordinateSystem
            Coordinate system of the image on which the gridded response is computed.
        time: double
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        freq: double
            Frequency in Hz
        station_index: int
            Station index, where 0 <= station_index < telescope.nr_stations
        field_index: int, optional
            Field index. Only relevant for dish telescopes such as VLA and ATCA. Default
            value is 0
            NOTE: field_index is a keyword-only argument

        Returns
        -------
        np.ndarray, np.complex64
            4d numpy array [height, width, 2, 2]
       )pbdoc",
          py::arg("coordinate_system"), py::arg("time"), py::arg("freq"),
          py::arg("station_index"), py::kw_only(), py::arg("field_index") = 0)
      .def(
          "gridded_response",
          [](Telescope &self, const CoordinateSystem &coordinate_system,
             double time, double freq,
             size_t field_idx) -> py::array_t<std::complex<float>> {
            const size_t nr_stations = self.GetNrStations();
            std::unique_ptr<GriddedResponse> grid_response =
                self.GetGriddedResponse(coordinate_system);
            std::vector<std::complex<float>> buffer(
                grid_response->GetStationBufferSize(nr_stations));
            grid_response->FullResponseAllStations(buffer.data(), time, freq,
                                                   field_idx);
            return cast_tensor<float>(buffer.data(),
                                      {nr_stations, coordinate_system.height,
                                       coordinate_system.width, 2, 2});
          },
          R"pbdoc(
        Compute the gridded response for all stations

        Parameters
        ----------
        coordinate_system: everybeam.CoordinateSystem
            Coordinate system of the image on which the gridded response is computed.
        time: double
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        freq: double
            Frequency in Hz
        field_index: int, optional
            Field index. Only relevant for dish telescopes such as VLA and ATCA. Default
            value is 0
            NOTE: field_index is a keyword-only argument

        Returns
        -------
        np.ndarray, np.complex64
            5d numpy array [nr_stations, height, width, 2, 2]
       )pbdoc",
          py::arg("coordinate_system"), py::arg("time"), py::arg("freq"),
          py::kw_only(), py::arg("field_index") = 0)
      .def(
          "undersampled_response",
          [](Telescope &self, const CoordinateSystem &coordinate_system,
             const py::array_t<double> &time, double freq,
             size_t undersampling_factor,
             const py::array_t<double> &baseline_weights, size_t field_id) {
            std::unique_ptr<GriddedResponse> grid_response =
                self.GetGriddedResponse(coordinate_system);
            std::vector<double> buffer(
                grid_response->GetIntegratedBufferSize());

            // Copy numpy to std::vector
            // NOTE: this could be optimized by making it a "view" rather than
            // a copy
            std::vector<double> bw_vec(baseline_weights.size());
            std::copy_n(baseline_weights.data(), baseline_weights.size(),
                        bw_vec.data());

            std::vector<double> time_vec(time.size());
            std::copy_n(time.data(), time_vec.size(), time_vec.data());

            // Call overload depending on the size of the time_vec
            if (time_vec.size() == 1) {
              grid_response->IntegratedFullResponse(
                  buffer.data(), time_vec[0], freq, field_id,
                  undersampling_factor, bw_vec);
            } else {
              grid_response->IntegratedFullResponse(
                  buffer.data(), time_vec, freq, field_id, undersampling_factor,
                  bw_vec);
            }
            // The Hermitian matrices need to be stored as a Matrix4x4
            const size_t npixels =
                coordinate_system.height * coordinate_system.width;

            return cast_tensor_mc<16, aocommon::MC4x4>(
                hmc_to_mc(buffer.data(), npixels),
                {coordinate_system.height, coordinate_system.height, 4, 4});
          },
          R"pbdoc(
        Compute the gridded response on an undersampled grid and (FFT) interpolate
        the result to the original grid.

        Parameters
        ----------
        coordinate_system: everybeam.CoordinateSystem
            Coordinate system of the image on which the gridded response is computed.
        time: double, np.1darray
            (Vector of) time(s) in modified Julian date, UTC, in seconds (MJD(UTC), s)
        freq: double
            Frequency in Hz
        undersampling_factor : int
            Undersampling factor, i.e. the coarsening factor between the original grid and the
            coarse resolution grid on which the beam will be evaluated.
        baseline_weights : np.1darray
            Vector containing the weights per baseline. Should have size (time.size() * nr_baselines),
            where nr_baselines equal telescope.nr_stations * (telescope.nr_stations + 1) // 2
        field_index: int, optional
            Field index. Only relevant for dish telescopes such as VLA and ATCA. Default
            value is 0
            NOTE: field_index is a keyword-only argument

        Returns
        -------
        np.ndarray, np.complex64
            4d numpy array [height, width, 4, 4], i.e. a Mueller matrix
            for every pixel in the image
       )pbdoc",
          py::arg("coordinate_system"), py::arg("time"), py::arg("freq"),
          py::arg("undersampling_factor"), py::arg("baseline_weights"),
          py::kw_only(), py::arg("field_index") = 0);

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
            return cast_tensor_mc<4, aocommon::MC2x2>(
                response, {nr_stations, nr_channels, 2, 2});
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
            return cast_tensor_mc<4, aocommon::MC2x2>(response,
                                                      {nr_channels, 2, 2});
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
            return cast_matrix(response);
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
            return cast_matrix(response);
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
             const py::array_t<double> &pydirection,
             const py::array_t<double> &pystation0,
             const py::array_t<double> &pytile0,
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

              return cast_matrix(response_diff_beam);
            } else {
              return cast_matrix(response);
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
              return cast_matrix(response_diff_beam);
            } else {
              return cast_matrix(response);
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
              return cast_matrix(response_diff_beam);
            } else {
              return cast_matrix(response);
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
              return cast_matrix(response_diff_beam);
            } else {
              return cast_matrix(response);
            }
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
            return cast_matrix(response);
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
            return cast_matrix(response);
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
