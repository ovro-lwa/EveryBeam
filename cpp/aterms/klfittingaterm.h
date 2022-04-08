// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_KL_FITTING_ATERM_H_
#define EVERYBEAM_ATERMS_KL_FITTING_ATERM_H_

#include "atermbase.h"
#include "../coords/coordutils.h"

namespace everybeam {
namespace aterms {

/**
 * Class that reads in H5Parm solution files and
 * fits Karhunen-Loève base functions to them\
 */
class KlFittingATerm final : public ATermBase {
 public:
  KlFittingATerm(const coords::CoordinateSystem& coordinate_system);

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
  bool Calculate(std::complex<float>* buffer, double time, double frequency,
                 size_t field_id, const double* uvw_in_m) override;

  double AverageUpdateTime() const override;

 private:
  coords::CoordinateSystem coordinate_system_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
