// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "fitsaterm.h"

using everybeam::aterms::FitsATerm;

FitsATerm::FitsATerm(size_t nAntenna,
                     const aocommon::CoordinateSystem& coordinate_system,
                     size_t max_support)
    : FitsATermBase(nAntenna, coordinate_system, max_support) {}

FitsATerm::~FitsATerm() = default;

void FitsATerm::OpenTECFiles(const std::vector<std::string>& filenames) {
  mode_ = Mode::kTEC;
  readers_.reserve(filenames.size());
  for (const std::string& filename : filenames) {
    readers_.emplace_back(filename, true, true);
    if (readers_.back().NFrequencies() != 1) {
      throw std::runtime_error(
          "FITS file for TEC A-terms has multiple frequencies in it");
    }
  }
  InitializeFromFiles(readers_);
}

void FitsATerm::OpenDiagGainFiles(const std::vector<std::string>& filenames) {
  mode_ = Mode::kDiagonal;
  readers_.reserve(filenames.size());
  for (const std::string& filename : filenames) {
    readers_.emplace_back(filename, true, true);
    if (readers_.back().NMatrixElements() != 4) {
      throw std::runtime_error(
          "FITS file for diagonal gains did not have 4 matrix elements in it");
    }
  }
  InitializeFromFiles(readers_);
}

bool FitsATerm::Calculate(std::complex<float>* buffer, double time,
                          double frequency, size_t, const double*) {
  size_t time_index;
  bool requires_recalculation;
  bool position_changed = FindFilePosition(buffer, time, frequency, time_index,
                                           requires_recalculation);
  if (!position_changed) {
    return false;
  } else {
    if (requires_recalculation) {
      ReadImages(buffer, time_index, frequency);
      StoreInCache(frequency, buffer);
    }
    return true;
  }
}

void FitsATerm::ReadImages(std::complex<float>* buffer, size_t time_index,
                           double frequency) {
  const size_t freq_index =
      round((frequency - readers_.front().FrequencyDimensionStart()) /
            readers_.front().FrequencyDimensionIncr());
  const size_t img_index =
      GetTimestep(time_index).img_index * NFrequencies() + freq_index;
  aocommon::FitsReader& reader = readers_[GetTimestep(time_index).reader_index];
  scratch_a_.resize(GetResampler().ScratchASize());
  scratch_b_.resize(GetResampler().ScratchBSize(reader));
  // TODO do this in parallel. Needs to fix Resampler too, as currently it can't
  // run in parallel when a window is used.
  for (size_t antenna_index = 0; antenna_index != NAntennas();
       ++antenna_index) {
    // In case there is only one antenna in the measurement set, copy it
    // to all antennas. This is not very efficient as the single image
    // is still read + resampled NAnt times, but it's a border case.
    const size_t antenna_file_index =
        img_index * NAntennas() +
        ((reader.NAntennas() == 1) ? 0 : antenna_index);
    std::complex<float>* antenna_buffer =
        buffer + antenna_index * Width() * Height() * 4;

    switch (mode_) {
      case Mode::kTEC: {
        // TODO When we are in the same timestep but at a different frequency,
        // it would be possible to skip reading and resampling, and immediately
        // call EvaluateTEC() with the "scratch" data still there.
        ReadAndResample(reader, antenna_file_index, scratch_a_, scratch_b_,
                        1.0);
        EvaluateTEC(antenna_buffer, scratch_b_.data(), frequency);
        break;
      }
      case Mode::kDiagonal: {
        size_t file_index = antenna_file_index * 4;
        for (size_t p = 0; p != 2; ++p) {
          ReadAndResample(reader, file_index, scratch_a_, scratch_b_, 1.0);
          CopyToRealPolarization(antenna_buffer, scratch_b_.data(), p * 3);

          ReadAndResample(reader, file_index + 1, scratch_a_, scratch_b_, 1.0);
          CopyToImaginaryPolarization(antenna_buffer, scratch_b_.data(), p * 3);

          file_index += 2;
        }
        SetPolarization(antenna_buffer, 1, std::complex<float>(0.0, 0.0));
        SetPolarization(antenna_buffer, 2, std::complex<float>(0.0, 0.0));
        break;
      }
    }
  }
}

void FitsATerm::EvaluateTEC(std::complex<float>* dest, const float* source,
                            double frequency) {
  for (size_t pixel = 0; pixel != Width() * Height(); ++pixel) {
    dest[pixel * 4] =
        std::polar(1.0, source[pixel] * -8.44797245e9 / frequency);
    dest[pixel * 4 + 1] = 0.0;
    dest[pixel * 4 + 2] = 0.0;
    dest[pixel * 4 + 3] = dest[pixel * 4];
  }
}

void FitsATerm::CopyToRealPolarization(std::complex<float>* dest,
                                       const float* source, size_t polIndex) {
  dest += polIndex;
  for (size_t i = 0; i != Width() * Height(); ++i) {
    dest[i * 4].real(source[i]);
  }
}

void FitsATerm::CopyToImaginaryPolarization(std::complex<float>* dest,
                                            const float* source,
                                            size_t polIndex) {
  dest += polIndex;
  for (size_t i = 0; i != Width() * Height(); ++i) {
    dest[i * 4].imag(source[i]);
  }
}

void FitsATerm::SetPolarization(std::complex<float>* dest, size_t polIndex,
                                std::complex<float> value) {
  dest += polIndex;
  for (size_t i = 0; i != Width() * Height(); ++i) {
    dest[i * 4] = value;
  }
}
