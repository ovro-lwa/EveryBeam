// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "griddedresponse.h"
#include "../common/fftresampler.h"

#include "../telescope/telescope.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/matrix4x4.h>
#include <vector>
#include <complex>
#include <stdexcept>

using aocommon::HMC4x4;
using aocommon::MC2x2;
using aocommon::MC4x4;
using aocommon::UVector;
using everybeam::common::FFTResampler;
using everybeam::griddedresponse::GriddedResponse;

void GriddedResponse::IntegratedResponse(
    BeamMode beam_mode, double* buffer, double time, double frequency,
    size_t field_id, size_t undersampling_factor,
    const std::vector<double>& baseline_weights) {
  const size_t nstations = telescope_->GetNrStations();
  const size_t nbaselines = nstations * (nstations + 1) / 2;
  if (baseline_weights.size() != nbaselines) {
    throw std::runtime_error("baseline_weights vector has incorrect size.");
  }
  // Scaling factor
  double baseline_total_weight =
      std::accumulate(baseline_weights.begin(), baseline_weights.end(), 0.0);

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
  size_t npixels = width_ * height_;
  std::vector<HMC4x4> matrices(npixels, HMC4x4::Zero());
  MakeIntegratedSnapshot(beam_mode, matrices, time, frequency, field_id,
                         baseline_weights.data());

  for (HMC4x4& matrix : matrices) {
    matrix /= baseline_total_weight;
  }

  DoFFTResampling(buffer, width_, height_, width_original, height_original,
                  matrices);

  // Reset coordinate members to original values
  width_ = width_original;
  height_ = height_original;
  dl_ = dl_original;
  dm_ = dm_original;
}

void GriddedResponse::IntegratedResponse(
    BeamMode beam_mode, double* buffer, const std::vector<double>& time_array,
    double frequency, size_t field_id, size_t undersampling_factor,
    const std::vector<double>& baseline_weights) {
  size_t nstations = telescope_->GetNrStations();
  size_t nbaselines = nstations * (nstations + 1) / 2;
  if (baseline_weights.size() != time_array.size() * nbaselines) {
    throw std::runtime_error("baseline_weights vector has incorrect size.");
  }

  // Scaling factor
  double baseline_total_weight =
      std::accumulate(baseline_weights.begin(), baseline_weights.end(), 0.0);

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

  for (std::size_t tstep = 0; tstep != time_array.size(); ++tstep) {
    MakeIntegratedSnapshot(beam_mode, matrices, time_array[tstep], frequency,
                           field_id,
                           baseline_weights.data() + tstep * nbaselines);
  }

  for (HMC4x4& matrix : matrices) {
    matrix /= baseline_total_weight;
  }

  DoFFTResampling(buffer, width_, height_, width_original, height_original,
                  matrices);

  // Reset coordinate members to original values
  width_ = width_original;
  height_ = height_original;
  dl_ = dl_original;
  dm_ = dm_original;
}

void GriddedResponse::MakeIntegratedSnapshot(
    BeamMode beam_mode, std::vector<HMC4x4>& matrices, double time,
    double frequency, size_t field_id,
    const double* baseline_weights_interval) {
  const size_t nstations = telescope_->GetNrStations();
  UVector<std::complex<float>> buffer_undersampled(
      GetStationBufferSize(nstations));
  ResponseAllStations(beam_mode, buffer_undersampled.data(), time, frequency,
                      field_id);

  const size_t npixels = width_ * height_;
  for (size_t y = 0; y != height_; ++y) {
    for (size_t x = 0; x != width_; ++x) {
      size_t index = 0;
      HMC4x4 gain = HMC4x4::Zero();
      for (size_t s1 = 0; s1 != nstations; ++s1) {
        size_t offset_s1 = (s1 * npixels + y * width_ + x) * 4;
        MC2x2 A(buffer_undersampled[offset_s1],
                buffer_undersampled[offset_s1 + 1],
                buffer_undersampled[offset_s1 + 2],
                buffer_undersampled[offset_s1 + 3]);
        for (size_t s2 = s1; s2 != nstations; ++s2) {
          size_t offset_s2 = offset_s1 + (s2 - s1) * npixels * 4;
          MC2x2 B(buffer_undersampled[offset_s2],
                  buffer_undersampled[offset_s2 + 1],
                  buffer_undersampled[offset_s2 + 2],
                  buffer_undersampled[offset_s2 + 3]);

          // Compute Mueller matrix and apply vec trick, see
          // https://en.wikipedia.org/wiki/Kronecker_product#Matrix_equations
          MC4x4 baseline_mueller =
              MC4x4::KroneckerProduct(B.HermTranspose().Transpose(), A) +
              MC4x4::KroneckerProduct(A.HermTranspose().Transpose(), B);

          // Note: keep scalar factor in brackets to avoid a redundant
          // element-wise multiplication of the HMC4x4 matrix. Factor 0.5 is to
          // account for Kronecker delta additions.
          gain += HMC4x4(baseline_mueller) *
                  (0.5 * baseline_weights_interval[index]);
          ++index;
        }
      }
      matrices[y * width_ + x] += gain;
    }
  }
}

void GriddedResponse::DoFFTResampling(
    double* buffer, int width_in, int height_in, int width_out, int height_out,
    const std::vector<aocommon::HMC4x4>& matrices) {
  // (FFT) resampling, run multi-threaded?
  FFTResampler resampler(width_in, height_in, width_out, height_out);
  resampler.SetWindowFunction(aocommon::WindowFunction::RaisedHann, true);
  UVector<float> lowres_input(width_in * height_in);
  // Loop over the "double"  representation of the HMC4x4 Hermitian matrix
  for (size_t p = 0; p != 16; ++p) {
    for (int i = 0; i != width_in * height_in; ++i) {
      lowres_input[i] = matrices[i].Data(p);
    }
    // Resample and write to the output buffer
    UVector<float> hires_output(width_out * height_out);
    resampler.Resample(lowres_input.data(), hires_output.data());
    // Copy hires_output (float) to buffer (double), conversion is implicit in
    // std::copy
    std::copy(hires_output.begin(), hires_output.end(),
              buffer + p * width_out * height_out);
  }
}
