// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_FOURIER_FITTER_
#define EVERYBEAM_ATERMS_FOURIER_FITTER_

#include <complex>
#include <utility>
#include <vector>

#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp>

namespace everybeam {
namespace aterms {

class [[gnu::visibility("default")]] FourierFitter {
 public:
  FourierFitter(std::size_t subgrid_size, std::size_t support,
                const std::vector<std::pair<float, float>>& directions);

  void Evaluate(const std::vector<std::complex<float>>& solutions,
                std::complex<float>* buffer) const;

 private:
  std::size_t subgrid_size_;
  std::size_t support_;
  xt::xtensor<std::complex<float>, 2> a_inv_;
};

}  // namespace aterms
}  // namespace everybeam
#endif
