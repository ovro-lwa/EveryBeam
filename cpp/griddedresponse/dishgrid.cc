// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "dishgrid.h"
#include "../telescope/dish.h"
#include "../circularsymmetric/voltagepattern.h"
#include "../circularsymmetric/vlabeam.h"

#include <aocommon/uvector.h>
#include <aocommon/matrix2x2.h>
#include <aocommon/matrix4x4.h>
#include <aocommon/hmatrix4x4.h>
#include <algorithm>

using aocommon::HMC4x4;
using aocommon::UVector;

namespace everybeam {
namespace griddedresponse {

void DishGrid::FullResponse(std::complex<float>* buffer, double,
                            double frequency,
                            [[maybe_unused]] size_t station_idx,
                            size_t field_id) {
  const telescope::Dish& dishtelescope =
      static_cast<const telescope::Dish&>(*telescope_);

  const double pdir_ra =
      dishtelescope.ms_properties_.field_pointing[field_id].first;
  const double pdir_dec =
      dishtelescope.ms_properties_.field_pointing[field_id].second;
  const std::array<double, 5> coefs =
      circularsymmetric::VLABeam::GetCoefficients("", frequency);
  const aocommon::UVector<double> coefs_vec(coefs.begin(), coefs.end());
  circularsymmetric::VoltagePattern vp(frequency, 53.0);
  vp.EvaluatePolynomial(coefs_vec, false);
  vp.Render(buffer, width_, height_, dl_, dm_, ra_, dec_, pdir_ra, pdir_dec,
            phase_centre_dl_, phase_centre_dm_, frequency);
}

void DishGrid::FullResponseAllStations(std::complex<float>* buffer, double,
                                       double frequency, size_t field_id) {
  FullResponse(buffer, 0.0, frequency, 0, field_id);

  const size_t station_buffer = width_ * height_ * 4;
  // Just repeat nstations times
  for (size_t i = 1; i != telescope_->GetNrStations(); ++i) {
    std::copy_n(buffer, station_buffer, buffer + i * station_buffer);
  }
}

void DishGrid::IntegratedFullResponse(double* buffer, double, double frequency,
                                      size_t field_id,
                                      size_t undersampling_factor,
                                      const std::vector<double>&) {
  // Copy coordinate members
  const size_t width_original = width_;
  const size_t height_original = height_;
  const double dl_original = dl_;
  const double dm_original = dm_;

  width_ /= undersampling_factor;
  height_ /= undersampling_factor;
  dl_ *= (double(width_original) / double(width_));
  dm_ *= (double(width_original) / double(width_));

  // Init (Hermitian) Mueller matrix for every pixel in the coarse grid
  const size_t npixels = width_ * height_;
  std::vector<HMC4x4> matrices(npixels, HMC4x4::Zero());
  MakeIntegratedDishSnapshot(matrices, frequency, field_id);

  DoFFTResampling(buffer, width_, height_, width_original, height_original,
                  matrices);

  // Reset coordinate members to original values
  width_ = width_original;
  height_ = height_original;
  dl_ = dl_original;
  dm_ = dm_original;
}

void DishGrid::MakeIntegratedDishSnapshot(
    std::vector<aocommon::HMC4x4>& matrices, double frequency,
    size_t field_id) {
  const size_t nstations = telescope_->GetNrStations();
  UVector<std::complex<float>> buffer_undersampled(
      GetStationBufferSize(nstations));
  // Assert that buffer size can accomodate Jones matrix on pixels
  assert(buffer_undersampled.size() >= width_ * height_ * 4);

  FullResponseAllStations(buffer_undersampled.data(), 0.0, frequency, field_id);
  // Loop over the pixels just once, and compute auto correlation
  for (size_t y = 0; y != height_; ++y) {
    for (size_t x = 0; x != width_; ++x) {
      const size_t offset = (y * width_ + x) * 4;
      const aocommon::MC2x2 A(&buffer_undersampled[offset]);
      // Compute Mueller matrix and apply vec trick, see
      // https://en.wikipedia.org/wiki/Kronecker_product#Matrix_equations
      // No need to add (+) and average (*0.5) for auto-correlation
      matrices[y * width_ + x] =
          HMC4x4::KroneckerProduct(A.HermTranspose().Transpose(), A);
    }
  }
}
}  // namespace griddedresponse
}  // namespace everybeam