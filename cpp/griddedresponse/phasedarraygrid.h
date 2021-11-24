// phasedarraygrid.h: Intermediate class for computing the gridded
// response for phased array radio telescopes such as LOFAR and SKA(-LOW)
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_PHASEDARRAYGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_PHASEDARRAYGRID_H_

#include "griddedresponse.h"
#include "../beamnormalisationmode.h"

namespace everybeam {
namespace griddedresponse {
class PhasedArrayGrid : public GriddedResponse {
 public:
  PhasedArrayGrid(const telescope::Telescope* telescope_ptr,
                  const coords::CoordinateSystem& coordinate_system);

  void Response(BeamMode beam_mode, std::complex<float>* buffer, double time,
                double frequency, size_t station_idx,
                size_t field_id) final override;

  void ResponseAllStations(BeamMode beam_mode, std::complex<float>* buffer,
                           double time, double frequency,
                           size_t field_id) override;

 protected:
  casacore::MDirection delay_dir_, tile_beam_dir_, preapplied_beam_dir_;
  vector3r_t station0_, tile0_, l_vector_itrf_, m_vector_itrf_, n_vector_itrf_,
      diff_beam_centre_;

  CorrectionMode preapplied_correction_mode_;
  BeamNormalisationMode beam_normalisation_mode_;
  bool use_channel_frequency_;
  double subband_frequency_;

 private:
  std::vector<aocommon::MC2x2F> inverse_central_gain_;

  std::vector<std::thread> threads_;

  struct Job {
    Job() {}
    Job(size_t y_, size_t antenna_idx_, size_t buffer_offset_)
        : y(y_), antenna_idx(antenna_idx_), buffer_offset(buffer_offset_) {}
    size_t y, antenna_idx, buffer_offset;
  };
  aocommon::Lane<Job>* lane_;

  bool CalculateBeamNormalisation(BeamMode beam_mode, double time,
                                  double frequency, size_t station_idx,
                                  aocommon::MC2x2F& inverse_gain) const;

  /**
   * @brief Method for computing the ITRF-vectors.
   * NOTE: method not thread-safe due to casacore dependencies.
   *
   * @param time
   */
  void SetITRFVectors(double time);

  void CalcThread(BeamMode beam_mode, bool apply_normalisation,
                  std::complex<float>* buffer, double time, double frequency);
};
}  // namespace griddedresponse
}  // namespace everybeam

#endif  // EVERYBEAM_GRIDDEDRESPONSE_PHASEDARRAYGRID_H_
