// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "beammode.h"
#include "load.h"
#include "station.h"
#include "griddedresponse/griddedresponse.h"
#include "pointresponse/phasedarraypoint.h"
#include "common/types.h"
#include "coords/coordutils.h"
#include "coords/itrfconverter.h"
#include "coords/itrfdirection.h"
#include "telescope/phasedarray.h"
#include "telescope/lofar.h"
#include "telescope/oskar.h"
#include "telescope/skamid.h"

#include <aocommon/matrix2x2.h>

namespace py = pybind11;

using casacore::MeasurementSet;

using aocommon::CoordinateSystem;
using everybeam::BeamMode;
using everybeam::BeamNormalisationMode;
using everybeam::Options;
using everybeam::vector3r_t;
using everybeam::coords::ITRFConverter;
using everybeam::coords::SetITRFVector;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::pointresponse::PhasedArrayPoint;
using everybeam::pointresponse::PointResponse;
using everybeam::telescope::LOFAR;
using everybeam::telescope::OSKAR;
using everybeam::telescope::PhasedArray;
using everybeam::telescope::SkaMid;
using everybeam::telescope::Telescope;

namespace {
// Convert pyarray of size 3 to vector3r_t
vector3r_t np2vector3r_t(const py::array_t<double>& pyarray) {
  auto r = pyarray.unchecked<1>();
  if (r.size() != 3) {
    throw std::runtime_error("Pyarray is of incorrect size, must be 3.");
  }
  return {r[0], r[1], r[2]};
}

// Cast complex buffer to numpy array
template <typename T>
py::array_t<std::complex<T>> cast_tensor(const std::complex<T>* buffer,
                                         const std::vector<size_t>& layout) {
  std::vector<ptrdiff_t> np_layout(layout.begin(), layout.end());
  return py::array_t<std::complex<T>>(np_layout, buffer);
}

// Cast aocommon::MC2x2 to py::array
py::array_t<std::complex<double>> cast_matrix(const aocommon::MC2x2& matrix) {
  // Solution from: https://github.com/pybind/pybind11/issues/1299
  return cast_tensor<double>(matrix.Data(), {2, 2});
}

// Cast vector of aocommon complex matrices (MC2x2/MC4x4) to numpy tensor
// with some additional size checks
template <size_t nelem, typename T>
py::array_t<std::complex<double>> cast_tensor_mc(
    const std::vector<T>& matrix, const std::vector<size_t>& layout) {
  size_t total_size = 1;
  for (size_t rank_size : layout) {
    total_size *= rank_size;
  }

  if (total_size != (matrix.size() * nelem)) {
    throw std::runtime_error("Casting mismatching shapes");
  }

  // Reinterpret_cast needed to flatten the nested std::array
  auto buffer = reinterpret_cast<const std::complex<double>*>(matrix.data());
  return cast_tensor<double>(buffer, layout);
}

// Convenience method throwing an error if provided station
// index exceeds number of stations
void check_station_index(size_t idx, size_t idx_max,
                         const std::string& prefix) {
  if (idx >= idx_max) {
    throw std::runtime_error(
        prefix + ": Requested station index exceeds number of stations.");
  }
}

// Convenience method to convert a buffer of [npixels * 16]
// (where 16 the number of entries of an aocmmon::HMC4x4 matrix
// to a vector of MC4x4 matrices
std::vector<aocommon::MC4x4> hmc_to_mc(float* buffer, size_t npixels) {
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
void apply_differential_beam(aocommon::MC2x2& lhs, const aocommon::MC2x2& rhs) {
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

  GriddedResponse* GetGriddedResponseWrapper(
      const CoordinateSystem& coordinate_system) {
    PYBIND11_OVERLOAD_PURE(GriddedResponse*, Telescope, GetGriddedResponse,
                           coordinate_system);
  }
};

/**
 * @brief Create a telescope object. As soon as the pybindings
 * support all the telescopes that are supported by the C++-interface,
 * this function becomes obsolete, and \c everybeam::Load() can be
 * used instead.
 *
 * @param name Path to MSet
 * @param options Options
 * @return std::unique_ptr<Telescope>
 */
std::unique_ptr<Telescope> create_telescope(const std::string& name,
                                            const everybeam::Options& options) {
  MeasurementSet ms(name);
  std::unique_ptr<Telescope> telescope;
  const everybeam::TelescopeType telescope_name =
      everybeam::GetTelescopeType(ms);
  switch (telescope_name) {
    case everybeam::TelescopeType::kAARTFAAC:
    case everybeam::TelescopeType::kLofarTelescope:
      telescope.reset(new LOFAR(ms, options));
      break;
    case everybeam::TelescopeType::kOSKARTelescope:
      telescope.reset(new OSKAR(ms, options));
      break;
    case everybeam::TelescopeType::kSkaMidTelescope:
      telescope.reset(new SkaMid(ms, options));
      break;
    default:
      throw std::runtime_error(
          "Currently, pybindings are only available for LOFAR and OSKAR "
          "MSets.");
  }
  return telescope;
}

std::unique_ptr<LOFAR> create_lofar(const std::string& name,
                                    const everybeam::Options& options) {
  return std::unique_ptr<LOFAR>{
      static_cast<LOFAR*>(create_telescope(name, options).release())};
}

std::unique_ptr<OSKAR> create_oskar(const std::string& name,
                                    const everybeam::Options& options) {
  return std::unique_ptr<OSKAR>{
      static_cast<OSKAR*>(create_telescope(name, options).release())};
}

std::unique_ptr<SkaMid> CreateSkaMid(const std::string& name,
                                     const everybeam::Options& options) {
  return std::unique_ptr<SkaMid>{
      static_cast<SkaMid*>(create_telescope(name, options).release())};
}

void init_telescope(py::module& m) {
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
        Number of stations in telescope.

        Returns
        -------
        int
       )pbdoc")
      .def_property_readonly("options", &Telescope::GetOptions,
                             R"pbdoc(
        Retrieve the specified options.

        Returns
        -------
        everybeam.Options
            Struct with options
       )pbdoc")
      // Documentation of overloaded function might not be parsed in the most
      // elegant way, see https://github.com/pybind/pybind11/issues/2619
      .def(
          "gridded_response",
          [](Telescope& self, const CoordinateSystem& coordinate_system,
             double time, double freq, size_t station_idx,
             size_t field_idx) -> py::array_t<std::complex<float>> {
            check_station_index(station_idx, self.GetNrStations(),
                                "gridded_response");
            py::array_t<std::complex<float>> buffer({coordinate_system.height,
                                                     coordinate_system.width,
                                                     size_t(2), size_t(2)});
            std::unique_ptr<GriddedResponse> grid_response =
                self.GetGriddedResponse(coordinate_system);
            grid_response->Response(BeamMode::kFull, buffer.mutable_data(),
                                    time, freq, station_idx, field_idx);
            return buffer;
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
            Station index, where ``0 <= station_index < telescope.nr_stations``
        field_index: int, optional
            Field index. Only relevant for dish telescopes such as VLA and ATCA. Default
            value is 0.

            **NOTE**: field_index is a keyword-only argument

        Returns
        -------
        np.ndarray, np.complex64
            4d numpy array ``[height, width, 2, 2]``
       )pbdoc",
          py::arg("coordinate_system"), py::arg("time"), py::arg("freq"),
          py::arg("station_index"), py::kw_only(), py::arg("field_index") = 0)
      .def(
          "gridded_response",
          [](Telescope& self, const CoordinateSystem& coordinate_system,
             double time, double freq,
             size_t field_idx) -> py::array_t<std::complex<float>> {
            const size_t nr_stations = self.GetNrStations();
            std::unique_ptr<GriddedResponse> grid_response =
                self.GetGriddedResponse(coordinate_system);
            py::array_t<std::complex<float>> buffer(
                {nr_stations, coordinate_system.height, coordinate_system.width,
                 size_t(2), size_t(2)});
            grid_response->ResponseAllStations(
                BeamMode::kFull, buffer.mutable_data(), time, freq, field_idx);
            return buffer;
          },
          R"pbdoc(
        Compute the gridded response for all stations.

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
            value is 0.

            **NOTE**: field_index is a keyword-only argument

        Returns
        -------
        np.ndarray, np.complex64
            5d numpy array ``[nr_stations, height, width, 2, 2]``
       )pbdoc",
          py::arg("coordinate_system"), py::arg("time"), py::arg("freq"),
          py::kw_only(), py::arg("field_index") = 0)
      .def(
          "undersampled_response",
          [](Telescope& self, const CoordinateSystem& coordinate_system,
             const py::array_t<double>& time, double freq,
             size_t undersampling_factor,
             const py::array_t<double>& baseline_weights, size_t field_id) {
            std::unique_ptr<GriddedResponse> grid_response =
                self.GetGriddedResponse(coordinate_system);
            std::vector<float> buffer(grid_response->GetIntegratedBufferSize());

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
              grid_response->IntegratedResponse(BeamMode::kFull, buffer.data(),
                                                time_vec[0], freq, field_id,
                                                undersampling_factor, bw_vec);
            } else {
              grid_response->IntegratedResponse(BeamMode::kFull, buffer.data(),
                                                time_vec, freq, field_id,
                                                undersampling_factor, bw_vec);
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
            value is 0.

            **NOTE**: field_index is a keyword-only argument

        Returns
        -------
        np.ndarray, np.complex64
            4d numpy array ``[height, width, 4, 4]``, i.e. a Mueller matrix
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
            Number of channels
       )pbdoc")
      .def(
          "station_name",
          [](PhasedArray& self, size_t idx) {
            if (idx >= self.GetNrStations()) {
              throw std::runtime_error(
                  "Requested station index exceeds number of stations.");
            }
            return self.GetStation(idx).GetName();
          },
          R"pbdoc(
        Get the station name given a station index.

        Parameters
        ----------
        station_index: int
            Station index

        Returns
        -------
        str
        )pbdoc",
          py::arg("station_index"))
      .def(
          "station_element_response",
          [](PhasedArray& self, size_t idx) {
            if (idx >= self.GetNrStations()) {
              throw std::runtime_error(
                  "Requested station index exceeds number of stations.");
            }
            return self.GetStation(idx).GetElementResponse();
          },
          R"pbdoc(
        Get the element response for a station, given a station index.

        Parameters
        ----------
        station_index: int
            Station index

        Returns
        -------
        An element response object for the station.
        )pbdoc",
          py::arg("station_index"))
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
      // Corresponds to evaluate3 in lofarbeam
      .def(
          "station_response",
          [](PhasedArray& self, double time, size_t idx, double freq,
             bool rotate) -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "station_response");
            vector3r_t direction;
            ITRFConverter itrf_converter(time);
            SetITRFVector(
                itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                direction);

            std::unique_ptr<PointResponse> point_response =
                self.GetPointResponse(time);
            PhasedArrayPoint& phased_array_point =
                static_cast<PhasedArrayPoint&>(*point_response);
            phased_array_point.SetParalacticRotation(rotate);

            std::mutex mutex;
            const aocommon::MC2x2 response =
                (self.GetOptions().beam_normalisation_mode ==
                 BeamNormalisationMode::kPreApplied)
                    ? aocommon::MC2x2::Unity()
                    : phased_array_point.Response(BeamMode::kFull, idx, freq,
                                                  direction, &mutex);
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
      // Corresponds to evaluate2 in lofarbeam
      .def(
          "station_response",
          [](PhasedArray& self, double time, size_t idx, size_t channel,
             bool rotate) -> py::array_t<std::complex<double>> {
            if (channel >= self.GetNrChannels()) {
              throw std::runtime_error(
                  "Requested channel index exceeds channel count.");
            }
            const double freq = self.GetChannelFrequency(channel);
            py::object py_station_response =
                py::cast(self).attr("station_response");
            return py_station_response(time, idx, freq, rotate);
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
            Apply paralactic angle rotation? ``[True/False]`` Defaults to ``True``

        Returns
        -------
        np.ndarray
            Response (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("channel_idx"),
          py::arg("rotate") = true)
      // Corresponds to evaluate1 in lofarbeam
      .def(
          "station_response",
          [](PhasedArray& self, double time, size_t idx,
             bool rotate) -> py::array_t<std::complex<double>> {
            // NOTE: reusing another station_response can be expected
            // to be slightly less efficient than making a dedicated
            // implementation - no use can be made of cached direction
            // vectors. It, however, considerably reduces the code base.
            const size_t nr_channels = self.GetNrChannels();

            py::array_t<std::complex<double>> response(
                {nr_channels, size_t(2), size_t(2)});
            py::object obj = py::cast(self);
            py::object py_station_response = obj.attr("station_response");
            for (size_t chan = 0; chan < nr_channels; ++chan) {
              response[py::make_tuple(chan, py::ellipsis())] =
                  py_station_response(time, idx, chan, rotate);
            }
            return response;
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
            Apply paralactic angle rotation? ``[True/False]`` Defaults to ``True``

        Returns
        -------
        np.ndarray
            rank 3 numpy array of shape ``[nr_channels, 2, 2]``
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("rotate") = true)
      // Corresponds to evaluate0 in lofarbeam
      .def(
          "station_response",
          [](PhasedArray& self, double time,
             bool rotate) -> py::array_t<std::complex<double>> {
            const size_t nr_stations = self.GetNrStations();
            const size_t nr_channels = self.GetNrChannels();

            py::object py_station_response =
                py::cast(self).attr("station_response");
            py::array_t<std::complex<double>> response(
                {nr_stations, nr_channels, size_t(2), size_t(2)});

            for (size_t sidx = 0; sidx < nr_stations; ++sidx) {
              for (size_t cidx = 0; cidx < nr_channels; ++cidx) {
                response[py::make_tuple(sidx, cidx, py::ellipsis())] =
                    py_station_response(time, sidx, cidx, rotate);
              }
            }
            return response;
          },
          R"pbdoc(
        Get station response in beam former direction for all stations and channels.

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        rotate: bool, optional
            Apply paralactic angle rotation? ``[True/False]`` Defaults to ``True``

        Returns
        -------
        np.ndarray
            rank 4 numpy array of shape ``[nr_stations, nr_channels, 2, 2]``
       )pbdoc",
          py::arg("time"), py::arg("rotate") = true)
      // Corresponds to evaluate4 from lofarbeam
      .def(
          "station_response",
          [](PhasedArray& self, double time, size_t idx, double freq,
             const py::array_t<double>& pydirection,
             const py::array_t<double>& pystation0,
             const py::array_t<double>& pytile0,
             bool rotate) -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "station_response");

            const vector3r_t direction = np2vector3r_t(pydirection);
            const vector3r_t station0 = np2vector3r_t(pystation0);
            const vector3r_t tile0 = np2vector3r_t(pytile0);

            std::unique_ptr<PointResponse> point_response =
                self.GetPointResponse(time);
            PhasedArrayPoint& phased_array_point =
                static_cast<PhasedArrayPoint&>(*point_response);
            phased_array_point.SetParalacticRotation(rotate);
            aocommon::MC2x2 response = phased_array_point.UnnormalisedResponse(
                BeamMode::kFull, idx, freq, direction, station0, tile0);

            // We can't delegate the normalisation to Response, in this case.
            // Since the diff_beam_centre direction would not be computed
            // correctly.
            if (self.GetOptions().beam_normalisation_mode ==
                BeamNormalisationMode::kPreApplied) {
              vector3r_t diff_beam_centre;
              ITRFConverter itrf_converter(time);
              SetITRFVector(
                  itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                  diff_beam_centre);

              aocommon::MC2x2 response_diff_beam =
                  phased_array_point.UnnormalisedResponse(
                      BeamMode::kFull, idx, freq, diff_beam_centre, station0,
                      tile0);
              apply_differential_beam(response_diff_beam, response);
              return cast_matrix(response_diff_beam);
            } else {
              return cast_matrix(response);
            }
            return cast_matrix(response);
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
            Apply paralactic angle rotation? ``[True/False]`` Defaults to ``True``

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
          [](PhasedArray& self, double time, size_t idx, double freq,
             const py::array_t<double> pydirection,
             const py::array_t<double> pystation0,
             bool rotate) -> py::array_t<std::complex<double>> {
            py::object py_station_response =
                py::cast(self).attr("station_response");
            return py_station_response(time, idx, freq, pydirection, pystation0,
                                       pystation0, rotate);
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
            Apply paralactic angle rotation? ``[True/False]`` Defaults to ``True``

        Returns
        -------
        np.ndarray
            Response (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("freq"),
          py::arg("direction"), py::arg("station0_direction"),
          py::arg("rotate") = true)
      .def(
          "station_response",
          [](PhasedArray& self, double time, size_t idx, double freq, double ra,
             double dec, size_t field_id) -> py::array_t<std::complex<float>> {
            std::unique_ptr<PointResponse> point_response =
                self.GetPointResponse(time);
            PhasedArrayPoint& phased_array_point =
                static_cast<PhasedArrayPoint&>(*point_response);

            py::array_t<std::complex<float>> response({size_t(2), size_t(2)});
            phased_array_point.Response(BeamMode::kFull,
                                        response.mutable_data(), ra, dec, freq,
                                        idx, field_id);
            return response;
          },
          R"pbdoc(
        Get station response in user-specified (ra, dec) direction. Delay direction is
        directly inferred from the provided MSet.

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        station_idx: int
            station index
        freq: float
            Frequency of the plane wave. (Hz)
        ra: float
            Right ascension coordinate of point of interest. (rad)
        dec: float
            Declination coordinate of point of interest. (rad)
        field_id: bool, optional
            Field index. (defaults to 0)

        Returns
        -------
        np.2darray
            Response (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("freq"),
          py::arg("ra"), py::arg("dec"), py::arg("field_id") = 0)
      .def(
          "element_response",
          [](PhasedArray& self, double time, size_t idx, size_t element_idx,
             double freq, const py::array_t<double> pydirection, bool is_local,
             bool rotate) -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "element_response");
            // TODO: need such a check for element index too
            const vector3r_t direction = np2vector3r_t(pydirection);

            std::unique_ptr<PointResponse> point_response =
                self.GetPointResponse(time);
            PhasedArrayPoint& phased_array_point =
                static_cast<PhasedArrayPoint&>(*point_response);
            phased_array_point.SetUseLocalCoordinateSystem(is_local);
            phased_array_point.SetParalacticRotation(rotate);

            const aocommon::MC2x2 response = phased_array_point.ElementResponse(
                idx, freq, direction, element_idx);

            if (self.GetOptions().beam_normalisation_mode ==
                BeamNormalisationMode::kPreApplied) {
              vector3r_t diff_beam_centre;
              ITRFConverter itrf_converter(time);
              SetITRFVector(
                  itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                  diff_beam_centre);
              aocommon::MC2x2 response_diff_beam =
                  phased_array_point.ElementResponse(
                      idx, freq, diff_beam_centre, element_idx);
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
            system is assumed. ``[True/False]`` Defaults to ``False``.
        rotate: bool, optional
            Apply paralactic angle rotation? ``[True/False]`` Defaults to ``True``

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
          [](PhasedArray& self, double time, size_t idx, double freq,
             const py::array_t<double> pydirection, bool is_local,
             bool rotate) -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "element_response");
            // TODO: need such a check for element index too
            const vector3r_t direction = np2vector3r_t(pydirection);

            std::unique_ptr<PointResponse> point_response =
                self.GetPointResponse(time);
            PhasedArrayPoint& phased_array_point =
                static_cast<PhasedArrayPoint&>(*point_response);
            phased_array_point.SetUseLocalCoordinateSystem(is_local);
            phased_array_point.SetParalacticRotation(rotate);

            // Avoid any beam normalisations, so compute station0 and tile0
            // manually
            ITRFConverter itrf_converter(time);
            vector3r_t station0;
            vector3r_t tile0;
            SetITRFVector(itrf_converter.ToDirection(self.GetDelayDirection()),
                          station0);
            SetITRFVector(
                itrf_converter.ToDirection(self.GetTileBeamDirection()), tile0);
            const aocommon::MC2x2 response =
                phased_array_point.UnnormalisedResponse(
                    BeamMode::kElement, idx, freq, direction, station0, tile0);

            if (self.GetOptions().beam_normalisation_mode ==
                BeamNormalisationMode::kPreApplied) {
              vector3r_t diff_beam_centre;
              SetITRFVector(
                  itrf_converter.ToDirection(self.GetPreappliedBeamDirection()),
                  diff_beam_centre);

              aocommon::MC2x2 response_diff_beam =
                  phased_array_point.UnnormalisedResponse(
                      BeamMode::kElement, idx, freq, diff_beam_centre, station0,
                      tile0);
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
            system is assumed. ``[True/False]`` Defaults to ``False``.
        rotate: bool, optional
            Apply paralactic angle rotation? ``[True/False]`` Defaults to ``True``

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
          [](PhasedArray& self, double time, size_t idx, double freq,
             const py::array_t<double> pydirection,
             const py::array_t<double> pystation0,
             const py::array_t<double> pytile0)
              -> py::array_t<std::complex<double>> {
            check_station_index(idx, self.GetNrStations(), "array_factor");
            const vector3r_t direction = np2vector3r_t(pydirection);
            const vector3r_t station0 = np2vector3r_t(pystation0);
            const vector3r_t tile0 = np2vector3r_t(pytile0);

            std::unique_ptr<PointResponse> point_response =
                self.GetPointResponse(time);
            PhasedArrayPoint& phased_array_point =
                static_cast<PhasedArrayPoint&>(*point_response);

            // Diagonal to 2x2 matrix
            const aocommon::MC2x2 response(
                phased_array_point.UnnormalisedResponse(BeamMode::kArrayFactor,
                                                        idx, freq, direction,
                                                        station0, tile0));
            return cast_matrix(response);
          },
          R"pbdoc(
        Get array factor for a given station in prescribed direction, with user-defined
        station0 and tile0 directions.

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
          py::arg("tile0_direction"))
      .def(
          "array_factor",
          [](PhasedArray& self, double time, size_t idx, double freq,
             const py::array_t<double> pydirection,
             const py::array_t<double> pystation0)
              -> py::array_t<std::complex<double>> {
            py::object py_array_factor = py::cast(self).attr("array_factor");
            return py_array_factor(time, idx, freq, pydirection, pystation0,
                                   pystation0);
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
          py::arg("direction"), py::arg("station0_direction"));

  py::class_<LOFAR, PhasedArray>(m, "LOFAR",
                                 R"pbdoc(
        Class to get beam responses for LOFAR observations.
        Inherits from :func:`~everybeam.PhasedArray`.
        )pbdoc")
      .def(py::init(&create_lofar),
           R"pbdoc(
        Initializes a LOFAR telescope.

        Parameters
        ----------
        ms: str
            Path to (LOFAR) Measurement Set
        options: everybeam.Options
            Struct specifying (beam) options for the provided
            Measurment Set
        )pbdoc",
           py::arg("ms"), py::arg("options"));

  py::class_<OSKAR, PhasedArray>(m, "OSKAR",
                                 R"pbdoc(
        Class to get beam responses for (simulated) SKA-LOW observations.
        Inherits from :func:`~everybeam.PhasedArray`.
        )pbdoc")
      .def(py::init(&create_oskar),
           R"pbdoc(
        Initializes an OSKAR telescope.

        Parameters
        ----------
        ms: str
            Path to (LOFAR) Measurement Set
        options: everybeam.Options
            Struct specifying (beam) options for the provided
            Measurment Set
        )pbdoc",
           py::arg("ms"), py::arg("options"));

  // TODO: other telescopes:
  // py::class_<MWA, Telescope>(m, "MWA");

  // py::class_<Dish, Telescope>(m, "Dish");

  // py::class_<VLA, Dish>(m, "MWA");

  py::class_<SkaMid, Telescope>(m, "SkaMid", R"pbdoc(
        Class to get beam responses for (simulated) SKA-MID observations.
        Inherits from :func:`~everybeam.Telescope`.
        )pbdoc")
      .def(py::init(&CreateSkaMid),
           R"pbdoc(
        Initializes a SKA-MID telescope.

        Parameters
        ----------
        ms: str
            Path to (LOFAR) Measurement Set
        options: everybeam.Options
            Struct specifying (beam) options for the provided
            Measurment Set
        )pbdoc",
           py::arg("ms"), py::arg("options"))
      .def_property_readonly("diameter", &SkaMid::GetDiameter,
                             R"pbdoc(
        Returns
        -------
        float
            Diameter of SKA-MID dish (m)
       )pbdoc")
      .def_property_readonly("blockage", &SkaMid::GetDiameter,
                             R"pbdoc(
        Returns
        -------
        float
            Blockage of SKA-MID dish due to receiver (m)
       )pbdoc")
      .def(
          // NOTE: duplicate of one of the overloads of
          // PhasedArray::station_response. Migrate to
          // Telescope base class?
          "station_response",
          [](SkaMid& self, double time, size_t idx, double freq, double ra,
             double dec, size_t field_id) -> py::array_t<std::complex<float>> {
            std::unique_ptr<PointResponse> point_response =
                self.GetPointResponse(time);

            py::array_t<std::complex<float>> response({size_t(2), size_t(2)});
            point_response->Response(BeamMode::kFull, response.mutable_data(),
                                     ra, dec, freq, idx, field_id);
            return response;
          },
          R"pbdoc(
        Get station response in user-specified (ra, dec) direction. Delay direction is
        directly inferred from the provided MSet.

        Parameters
        ----------
        time: double
            Evaluation response at time.
            Time in modified Julian date, UTC, in seconds (MJD(UTC), s)
        station_idx: int
            station index
        freq: float
            Frequency of the plane wave. (Hz)
        ra: float
            Right ascension coordinate of point of interest. (rad)
        dec: float
            Declination coordinate of point of interest. (rad)
        field_id: bool, optional
            Field index. (defaults to 0)

        Returns
        -------
        np.2darray
            Response (Jones) matrix
       )pbdoc",
          py::arg("time"), py::arg("station_idx"), py::arg("freq"),
          py::arg("ra"), py::arg("dec"), py::arg("field_id") = 0);
}
