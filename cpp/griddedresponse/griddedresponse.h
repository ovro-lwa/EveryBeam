// griddedresponse.h: Base class for computing the (gridded) telescope response.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_GRIDDEDRESPONSE_H_
#define EVERYBEAM_GRIDDEDRESPONSE_GRIDDEDRESPONSE_H_

#include "./../coords/coordutils.h"
#include "./../coords/itrfdirection.h"
#include "./../coords/itrfconverter.h"
#include "../beammode.h"

#include <aocommon/hmatrix4x4.h>

#include <memory>
#include <vector>
#include <thread>

#include <casacore/measures/Measures/MDirection.h>

namespace everybeam {

namespace telescope {
class Telescope;
}

namespace griddedresponse {

/**
 * @brief Virtual base class to compute the gridded response
 *
 */
class GriddedResponse {
 public:
  virtual ~GriddedResponse() {}

  /**
   * @brief Compute the beam for a single station, given a prescribed beam mode.
   * Result is stored in the output buffer, which should accommodate a Jones
   * matrix (4 complex floats per pixel.
   *
   * @param beam_mode Selects beam mode (BeamMode::kElement,
   * BeamMode::kArrayFactor or BeamMode::kFull)
   * @param buffer Output buffer, compute and set size with
   * GriddedResponse::GetStationBufferSize(1)
   * @param station_idx Station index, must be smaller than number of stations
   * in the Telescope
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param frequency Frequency (Hz)
   */
  virtual void Response(BeamMode beam_mode, std::complex<float>* buffer,
                        double time, double freq, size_t station_idx,
                        size_t field_id) = 0;

  /**
   * @brief Compute the array factor for all stations in a Telescope.
   * Result is stored in the output buffer, which should accommodate a Jones
   * matrix (4 complex valued floats) per pixel for each station.
   *
   * @param beam_mode Selects beam mode (element, array factor or full)
   * @param buffer Output buffer, compute and set size with
   * GriddedResponse::GetStationBufferSize()
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param frequency Frequency (Hz)
   */
  virtual void ResponseAllStations(BeamMode beam_mode,
                                   std::complex<float>* buffer, double time,
                                   double frequency, size_t field_id) = 0;

  /**
   * @brief Calculate integrated/undersampled beam for a single time step.
   * This function makes use of @ref MakeIntegratedSnapshot(). Subclasses
   * may override MakeIntegratedSnapshot() to implement a more efficient
   * version.
   *
   * @param buffer Buffer for storing the result, should have size width *
   * height * 16
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param frequency Frequency (Hz)
   * @param field_id Field id
   * @param undersampling_factor Undersampling factor
   * @param baseline_weights Baseline weights, size should equal
   * Telescope::GetNrStations() *  (Telescope::GetNrStations() + 1)/2
   */

  /**
   * @brief Calculate integrated/undersampled beam for a single time step.
   * This function makes use of @ref MakeIntegratedSnapshot(). Subclasses
   * may override MakeIntegratedSnapshot() to implement a more efficient
   * version.
   *
   * @param beam_mode Selects beam mode (element, array factor or full)
   * @param buffer Buffer for storing the result, should have size width *
   * height * 16
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param frequency Frequency (Hz)
   * @param field_id Field id
   * @param undersampling_factor Undersampling factor
   * @param baseline_weights Baseline weights, size should equal
   * Telescope::GetNrStations() *  (Telescope::GetNrStations() + 1)/2
   */
  virtual void IntegratedResponse(BeamMode beam_mode, float* buffer,
                                  double time, double frequency,
                                  size_t field_id, size_t undersampling_factor,
                                  const std::vector<double>& baseline_weights);

  /**
   * @brief Calculate integrated beam over multiple time steps.
   * This function makes use of @ref MakeIntegratedSnapshot(). Subclasses
   * may override MakeIntegratedSnapshot() to implement a more efficient
   * version.
   *
   * @param beam_mode Selects beam mode (element, array factor or full)
   * @param buffer Buffer for storing the result, should have size width *
   * height * 16
   * @param time_array Vector with probing times, modified Julian date, UTC, in
   * seconds (MJD(UTC), s).
   * @param frequency Frequency (Hz)
   * @param field_id Field id
   * @param undersampling_factor Undersampling factor
   * @param baseline_weights Baseline weights, size should equal
   * (Telescope::GetNrStations() *  (Telescope::GetNrStations() + 1)/2) *
   * time_array.size()
   */
  virtual void IntegratedResponse(BeamMode beam_mode, float* buffer,
                                  const std::vector<double>& time_array,
                                  double frequency, size_t field_id,
                                  size_t undersampling_factor,
                                  const std::vector<double>& baseline_weights);

  std::size_t GetStationBufferSize(std::size_t nstations) const {
    return nstations * width_ * height_ * 4u;
  }

  std::size_t GetIntegratedBufferSize() const { return width_ * height_ * 16u; }

  /**
   * Allow undersampled calculation of the beam grid? If true, the beam
   * will be calculated on a smaller grid and FFT resampled to a larger
   * grid. This is important for arrays that produce large images and for
   * which the beam calculations are expensive, such as MWA, LOFAR and
   * SKA. On the other hand, beams that are fast to calculate and might
   * have sharp cut-offs, should not be undersampled.
   */
  virtual bool PerformUndersampling() const { return true; }

 protected:
  /**
   * @brief Construct a new Gridded Response object
   *
   * @param telescope_ptr Pointer to telescope::Telescope object
   * @param coordinate_system CoordinateSystem struct
   */
  GriddedResponse(const telescope::Telescope* telescope_ptr,
                  const coords::CoordinateSystem& coordinate_system)
      : telescope_(telescope_ptr),
        width_(coordinate_system.width),
        height_(coordinate_system.height),
        ra_(coordinate_system.ra),
        dec_(coordinate_system.dec),
        dl_(coordinate_system.dl),
        dm_(coordinate_system.dm),
        phase_centre_dl_(coordinate_system.phase_centre_dl),
        phase_centre_dm_(coordinate_system.phase_centre_dm){};

  static void DoFFTResampling(float* buffer, int width_in, int height_in,
                              int width_out, int height_out,
                              const std::vector<aocommon::HMC4x4>& matrices);

  const telescope::Telescope* telescope_;
  size_t width_, height_;
  double ra_, dec_, dl_, dm_, phase_centre_dl_, phase_centre_dm_;

 private:
  /**
   * @brief Calculate a baseline-integrated snapshot.
   * By default, this function will request the response for all antennas and
   * perform a weighted average. (Partly) homogenous arrays may implement a
   * faster implementation by overriding this method.
   */
  virtual void MakeIntegratedSnapshot(BeamMode beam_mode,
                                      std::vector<aocommon::HMC4x4>& matrices,
                                      double time, double frequency,
                                      size_t field_id,
                                      const double* baseline_weights_interval);
};
}  // namespace griddedresponse
}  // namespace everybeam
#endif  // EVERYBEAM_GRIDDEDRESPONSE_GRIDDEDRESPONSE_H_
