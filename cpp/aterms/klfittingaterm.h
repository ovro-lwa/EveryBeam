// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_KL_FITTING_ATERM_H_
#define EVERYBEAM_ATERMS_KL_FITTING_ATERM_H_

#include "atermbase.h"

#include <memory>

#include <schaapcommon/h5parm/h5parm.h>
#include <aocommon/coordinatesystem.h>

namespace everybeam {
namespace aterms {

class KlFitter;

/**
 * Class that reads in H5Parm solution files and
 * fits Karhunen-Loève base functions to them\
 */
class [[gnu::visibility("default")]] KlFittingATerm final : public ATermBase {
 public:
  KlFittingATerm(const std::vector<std::string>& station_names_ms,
                 const aocommon::CoordinateSystem& coordinate_system, int order,
                 bool use_phasor_fit = true);

  // The destructor is declared here and defined in klfittingaterm.cc, to
  // be the default destructor. This split is needed, because here the size of
  // forward declared class KlFitter is unknown.
  ~KlFittingATerm();

  /**
   * @brief Read h5parm file given a path. Subsequent Calculate() calls will use
   * the parameters in the h5parm file for aterm calculation
   *
   * @param filename path to h5parm file
   */
  void Open(const std::string& filename);

  /**
   * @brief
   *
   * @param buffer Buffer
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s)
   * @param frequency Freq (Hz)
   * @param field_id
   * @param uvw_in_m
   * @return true Results are updated
   * @return false No need to update the result, cached result can be used
   */
  bool Calculate(std::complex<float> * buffer, double time, double frequency,
                 size_t field_id, const double* uvw_in_m) override;

  /**
   * @brief Set the update interval
   *
   * @param update_interval Update interval (in s)
   */
  void SetUpdateInterval(double update_interval) {
    update_interval_ = update_interval;
  }

  /**
   * @brief Get average update time, fixed value as set by SetUpdateInterval()
   *
   * @return double
   */
  double AverageUpdateTime() const final override { return update_interval_; }

 private:
  std::vector<std::string> station_names_ms_;
  const aocommon::CoordinateSystem coordinate_system_;
  int order_;
  schaapcommon::h5parm::SolTab phase_soltab_;
  double update_interval_;
  std::unique_ptr<KlFitter> kl_fitter_;
  double last_aterm_update_;
  std::size_t nr_directions_;
  bool use_phasor_fit_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
