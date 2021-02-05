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
using aocommon::MC2x2;
using aocommon::MC4x4;
using aocommon::UVector;
using everybeam::griddedresponse::DishGrid;

void DishGrid::CalculateStation(std::complex<float>* buffer, double,
                                double frequency, size_t, size_t field_id) {
  const telescope::Dish& dishtelescope =
      static_cast<const telescope::Dish&>(*telescope_);

  double pdir_ra = dishtelescope.ms_properties_.field_pointing[field_id].first,
         pdir_dec =
             dishtelescope.ms_properties_.field_pointing[field_id].second;
  std::array<double, 5> coefs =
      circularsymmetric::VLABeam::GetCoefficients("", frequency);
  circularsymmetric::VoltagePattern vp(frequency, 53.0);
  aocommon::UVector<double> coefs_vec(coefs.begin(), coefs.end());
  vp.EvaluatePolynomial(coefs_vec, false);
  vp.Render(buffer, width_, height_, dl_, dm_, ra_, dec_, pdir_ra, pdir_dec,
            phase_centre_dl_, phase_centre_dm_, frequency);
}

void DishGrid::CalculateAllStations(std::complex<float>* buffer, double,
                                    double frequency, size_t field_id) {
  CalculateStation(buffer, 0., frequency, 0, field_id);

  // Just repeat nstations times
  for (size_t i = 1; i != telescope_->GetNrStations(); ++i) {
    std::copy_n(buffer, width_ * height_ * 4,
                buffer + i * width_ * height_ * 4);
  }
}

void DishGrid::CalculateIntegratedResponse(double* buffer, double,
                                           double frequency, size_t field_id,
                                           size_t undersampling_factor,
                                           const std::vector<double>&) {
  // Copy coordinate members
  size_t width_original = width_, height_original = height_;
  double dl_original = dl_, dm_original = dm_;

  width_ /= undersampling_factor;
  height_ /= undersampling_factor;
  dl_ *= (double(width_original) / double(width_));
  dm_ *= (double(width_original) / double(width_));

  // Init (Hermitian) Mueller matrix for every pixel in the coarse grid
  size_t npixels = width_ * height_;
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
  size_t nstations = telescope_->GetNrStations();
  UVector<std::complex<float>> buffer_undersampled(
      GetStationBufferSize(nstations));
  // Assert that buffer size can accomodate Jones matrix on pixels
  assert(buffer_undersampled.size() >= width_ * height_ * 4);

  CalculateAllStations(buffer_undersampled.data(), 0., frequency, field_id);
  // Loop over the pixels just once, and compute auto correlation
  for (size_t y = 0; y != height_; ++y) {
    for (size_t x = 0; x != width_; ++x) {
      size_t offset = (y * width_ + x) * 4;
      MC2x2 A(buffer_undersampled[offset], buffer_undersampled[offset + 1],
              buffer_undersampled[offset + 2], buffer_undersampled[offset + 3]);
      // Compute Mueller matrix and apply vec trick, see
      // https://en.wikipedia.org/wiki/Kronecker_product#Matrix_equations
      // No need to add (+) and average (*0.5) for auto-correlation
      matrices[y * width_ + x] =
          HMC4x4::KroneckerProduct(A.HermTranspose().Transpose(), A);
    }
  }
}
