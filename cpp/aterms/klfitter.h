// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_KL_FITTER_
#define EVERYBEAM_ATERMS_KL_FITTER_

#include <complex>
#include <utility>
#include <vector>

#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp>

namespace everybeam {
namespace aterms {

class [[gnu::visibility("default")]] KlFitter {
 public:
  KlFitter(std::size_t subgrid_size, int order,
           const std::vector<std::pair<float, float>>& directions);

  void Evaluate(const std::vector<float>& solutions, float* buffer) const;

 private:
  std::size_t subgrid_size_;
  xt::xtensor<float, 3> fitting_matrix_;
};

}  // namespace aterms
}  // namespace everybeam
#endif
