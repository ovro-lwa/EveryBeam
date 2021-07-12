// phasedarraygrid.h: Intermediate class for computing the gridded
// response for phased array radio telescopes such as LOFAR and SKA(-LOW)
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_PHASEDARRAYGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_PHASEDARRAYGRID_H_

#include "griddedresponse.h"

#include "../correctionmode.h"

namespace everybeam {
namespace griddedresponse {
class PhasedArrayGrid : public GriddedResponse {
 public:
  PhasedArrayGrid(const telescope::Telescope* telescope_ptr,
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
                        size_t field_id) final override;

  /**
   * @brief Compute the Beam response for all stations in a Telescope
   *
   * @param buffer Output buffer, compute and set size with
   * GriddedResponse::GetStationBufferSize()
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency (Hz)
   */
  void CalculateAllStations(std::complex<float>* buffer, double time,
                            double frequency, size_t field_id) final override;

 protected:
  casacore::MDirection delay_dir_, tile_beam_dir_, preapplied_beam_dir_;
  vector3r_t station0_, tile0_, l_vector_itrf_, m_vector_itrf_, n_vector_itrf_,
      diff_beam_centre_;

  bool use_differential_beam_;
  CorrectionMode preapplied_correction_mode_;
  bool use_channel_frequency_;
  double subband_frequency_;

 private:
  std::vector<aocommon::MC2x2F> inverse_central_gain_;

  std::size_t nthreads_;
  std::vector<std::thread> threads_;

  struct Job {
    Job() {}
    Job(size_t y_, size_t antenna_idx_, size_t buffer_offset_)
        : y(y_), antenna_idx(antenna_idx_), buffer_offset(buffer_offset_) {}
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
