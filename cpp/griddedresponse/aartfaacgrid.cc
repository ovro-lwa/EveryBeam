// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "aartfaacgrid.h"

namespace everybeam {
namespace griddedresponse {

void AartfaacGrid::ResponseAllStations(BeamMode beam_mode,
                                       std::complex<float>* buffer, double time,
                                       double frequency, size_t field_id) {
  Response(beam_mode, buffer, time, frequency, 0, field_id);
  const size_t station_buffer = width_ * height_ * 4;
  // Repeated copy for n_stations
  for (size_t i = 1; i != telescope_->GetNrStations(); ++i) {
    std::copy_n(buffer, station_buffer, buffer + i * station_buffer);
  }
}

void AartfaacGrid::MakeIntegratedSnapshot(
    BeamMode beam_mode, std::vector<aocommon::HMC4x4>& matrices, double time,
    double frequency, size_t field_id,
    const double* baseline_weights_interval) {
  const size_t n_stations = telescope_->GetNrStations();
  aocommon::UVector<std::complex<float>> buffer_undersampled(
      GetStationBufferSize(n_stations));
  ResponseAllStations(beam_mode, buffer_undersampled.data(), time, frequency,
                      field_id);

  // For Aartfaac, we can simply weight a (time) snapshot with the accumulated
  // baseline weights
  const size_t n_baselines = n_stations * (n_stations + 1) / 2;
  double snapshot_weight = 0.;
  for (size_t index = 0; index != n_baselines; ++index) {
    snapshot_weight += baseline_weights_interval[index];
  }

  for (size_t y = 0; y != height_; ++y) {
    for (size_t x = 0; x != width_; ++x) {
      size_t offset = (y * width_ + x) * 4;
      const aocommon::MC2x2 A(&buffer_undersampled[offset]);

      // Mueller matrix constant for all baselines, so just compute once for
      // each individual pixel
      matrices[y * width_ + x] =
          aocommon::HMC4x4::KroneckerProduct(A.HermTranspose().Transpose(), A) *
          snapshot_weight;
    }
  }
}
}  // namespace griddedresponse
}  // namespace everybeam
