// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_FOURIER_FITTING_ATERM_H_
#define EVERYBEAM_ATERMS_FOURIER_FITTING_ATERM_H_

#include <complex>
#include <map>
#include <memory>
#include <vector>

#include <aocommon/coordinatesystem.h>
#include <aocommon/uvector.h>
#include <schaapcommon/h5parm/h5parm.h>
#include <numeric>

#include "atermbase.h"
#include "cache.h"

namespace everybeam {
namespace aterms {

class FourierFitter;

/**
 * Class that reads in H5Parm solution files and
 * fits Fourier base functions to them,
 * leading to aterms that have both limited support and
 * an exact reconstruction at the points of interest.
 */
class [[gnu::visibility("default")]] FourierFittingATerm final
    : public ATermBase {
 public:
  /**
   * @brief Constructs FourierFittingATerm
   *
   * @param station_names_ms vector of station names to compute aterms for
   * @param coordinate_system coordinate system of the aterms
   * @param support support of the aterm in the Fourier domain
   */
  FourierFittingATerm(const std::vector<std::string>& station_names_ms,
                      const aocommon::CoordinateSystem& coordinate_system,
                      int support);

  // The destructor is declared here and defined in fourierfittingaterm.cc, to
  // be the default destructor. This split is needed, because here the size of
  // forward declared class FourierFitter is unknown.
  ~FourierFittingATerm();

  /**
   * @brief Read h5parm file given a path. Subsequent Calculate() calls will use
   * the parameters in the h5parm file for aterm calculation
   *
   * @param filename path to h5parm file
   */
  void Open(const std::string& filename);

  /**
   * @brief Calculate aterms and write the result into buffer
   *
   * @param buffer Buffer
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s)
   * @param frequency Freq (Hz) - not used at the moment
   * @param field_id Irrelevant for fourierfitting aterms
   * @param uvw_in_m Irrelevant for fourierfitting aterms
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
  std::unique_ptr<FourierFitter> fourier_fitter_;
  const std::vector<std::string> station_names_ms_;
  const aocommon::CoordinateSystem coordinate_system_;

  schaapcommon::h5parm::SolTab phase_soltab_;
  int support_;
  std::size_t nr_directions_;
  double update_interval_;
  double last_aterm_update_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
