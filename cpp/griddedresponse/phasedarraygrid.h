// phasedarraygrid.h: Intermediate class for computing the gridded
// response for phased array radio telescopes such as LOFAR and SKA(-LOW)
//
// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the EveryBeam software suite.
// The EveryBeam software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The EveryBeam software suite is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the EveryBeam software suite. If not, see
// <http://www.gnu.org/licenses/>.
//
// $Id$

#ifndef EVERYBEAM_GRIDDEDRESPONSE_PHASEDARRAYGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_PHASEDARRAYGRID_H_

#include "griddedresponse.h"

namespace everybeam {
namespace griddedresponse {
class PhasedArrayGrid : public GriddedResponse {
 public:
  PhasedArrayGrid(telescope::Telescope* telescope_ptr,
                  const coords::CoordinateSystem& coordinate_system);

  /**
   * @brief Compute the Beam response for a single station
   *
   * @param buffer Output buffer, compute and set size with
   * GriddedResponse::GetBufferSize(1)
   * @param station_idx Station index, must be smaller than number of stations
   * in the Telescope
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency (Hz)
   */
  void CalculateStation(std::complex<float>* buffer, double time,
                        double frequency, size_t station_idx,
                        size_t field_id) override;

  /**
   * @brief Compute the Beam response for all stations in a Telescope
   *
   * @param buffer Output buffer, compute and set size with
   * GriddedResponse::GetStationBufferSize()
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency (Hz)
   */
  void CalculateAllStations(std::complex<float>* buffer, double time,
                            double frequency, size_t field_id) override;

 protected:
  casacore::MDirection delay_dir_, tile_beam_dir_, preapplied_beam_dir_;
  vector3r_t station0_, tile0_, l_vector_itrf_, m_vector_itrf_, n_vector_itrf_,
      diff_beam_centre_;

  bool use_differential_beam_, use_channel_frequency_;
  double subband_frequency_;

 private:
  std::vector<aocommon::MC2x2F> inverse_central_gain_;

  std::size_t nthreads_;
  std::vector<std::thread> threads_;

  struct Job {
    size_t y, antenna_idx, buffer_offset;
  };
  aocommon::Lane<Job>* lane_;

  /**
   * @brief Method for computing the ITRF-vectors.
   * NOTE: method not thread-safe due to casacore dependencies.
   *
   * @param time
   */
  void SetITRFVectors(double time);

  void CalcThread(std::complex<float>* buffer, double time, double frequency);
};
}  // namespace griddedresponse
}  // namespace everybeam

#endif  // EVERYBEAM_GRIDDEDRESPONSE_PHASEDARRAYGRID_H_