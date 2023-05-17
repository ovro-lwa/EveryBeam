// phasedarraygrid.h: Intermediate class for computing the gridded
// response for phased array radio telescopes such as LOFAR and SKA(-LOW)
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_GRIDDEDRESPONSE_PHASEDARRAYGRID_H_
#define EVERYBEAM_GRIDDEDRESPONSE_PHASEDARRAYGRID_H_

#include "griddedresponse.h"
#include "../beamnormalisationmode.h"
#include "../phasedarrayresponse.h"

namespace aocommon {
template <typename Tp>
class Lane;
}

namespace everybeam {
namespace griddedresponse {
class [[gnu::visibility("default")]] PhasedArrayGrid
    : public GriddedResponse,
      protected PhasedArrayResponse {
 public:
  PhasedArrayGrid(const telescope::Telescope* telescope_ptr,
                  const aocommon::CoordinateSystem& coordinate_system);

  void Response(BeamMode beam_mode, std::complex<float> * buffer, double time,
                double frequency, size_t station_idx, size_t field_id)
      final override;

  void ResponseAllStations(BeamMode beam_mode, std::complex<float> * buffer,
                           double time, double frequency, size_t field_id)
      override;

 private:
  vector3r_t l_vector_itrf_;
  vector3r_t m_vector_itrf_;
  vector3r_t n_vector_itrf_;
  std::vector<aocommon::MC2x2F> inverse_central_gain_;

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

  void CalcThread(BeamMode beam_mode, bool apply_normalisation,
                  std::complex<float>* buffer, double time, double frequency);
};
}  // namespace griddedresponse
}  // namespace everybeam

#endif  // EVERYBEAM_GRIDDEDRESPONSE_PHASEDARRAYGRID_H_
