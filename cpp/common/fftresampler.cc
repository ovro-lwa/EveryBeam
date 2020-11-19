// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "fftresampler.h"

#include <complex>
#include <iostream>

using aocommon::WindowFunction;
using everybeam::common::FFTResampler;

FFTResampler::FFTResampler(size_t width_in, size_t height_in, size_t width_out,
                           size_t height_out)
    : width_in_(width_in),
      height_in_(height_in),
      width_out_(width_out),
      height_out_(height_out),
      fft_width_(std::max(width_in, width_out)),
      fft_height_(std::max(height_in, height_out)),
      window_function_(WindowFunction::Rectangular),
      tukey_inset_size_(0.0),
      correct_window_(false) {
  float* input_data = reinterpret_cast<float*>(
      fftwf_malloc(fft_width_ * fft_height_ * sizeof(float)));
  fftwf_complex* fft_data = reinterpret_cast<fftwf_complex*>(
      fftwf_malloc(fft_width_ * fft_height_ * sizeof(fftwf_complex)));
  in_to_f_plan_ = fftwf_plan_dft_r2c_2d(height_in_, width_in_, input_data,
                                        fft_data, FFTW_ESTIMATE);
  f_to_out_plan_ = fftwf_plan_dft_c2r_2d(height_out_, width_out_, fft_data,
                                         input_data, FFTW_ESTIMATE);
  fftwf_free(fft_data);
  fftwf_free(input_data);
}

FFTResampler::~FFTResampler() {
  fftwf_destroy_plan(in_to_f_plan_);
  fftwf_destroy_plan(f_to_out_plan_);
}

void FFTResampler::RunSingle(const Task& task, bool skip_window) const {
  float* ptr_end = task.input + width_in_ * height_in_;
  for (float* i = task.input; i != ptr_end; ++i) {
    if (!std::isfinite(*i)) *i = 0.0;
  }

  if (window_function_ != WindowFunction::Rectangular && !skip_window)
    ApplyWindow(task.input);

  size_t fft_in_width = width_in_ / 2 + 1;
  std::complex<float>* fft_data = reinterpret_cast<std::complex<float>*>(
      fftwf_malloc(fft_in_width * height_in_ * sizeof(std::complex<float>)));
  fftwf_execute_dft_r2c(in_to_f_plan_, task.input,
                        reinterpret_cast<fftwf_complex*>(fft_data));

  size_t fft_out_width = width_out_ / 2 + 1;
  // TODO this can be done without allocating more mem!
  std::complex<float>* fft_data_new = reinterpret_cast<std::complex<float>*>(
      fftwf_malloc(fft_out_width * height_out_ * sizeof(std::complex<float>)));
  std::uninitialized_fill_n(fft_data_new, fft_out_width * height_out_,
                            std::complex<float>(0));

  size_t old_mid_x = width_in_ / 2;
  size_t new_mid_x = width_out_ / 2;

  size_t min_width = std::min(width_in_, width_out_);
  size_t min_height = std::min(height_in_, height_out_);

  size_t min_mid_x = min_width / 2;
  size_t min_mid_y = min_height / 2;

  float factor = 1.0 / (min_width * min_height);

  for (size_t y = 0; y != min_height; ++y) {
    size_t y_old = y - min_mid_y + height_in_;
    size_t y_new = y - min_mid_y + height_out_;
    if (y_old >= height_in_) y_old -= height_in_;
    if (y_new >= height_out_) y_new -= height_out_;

    // The last dimension is stored half
    for (size_t x = 0; x != min_mid_x; ++x) {
      size_t index_old = x + y_old * (old_mid_x + 1);
      size_t index_new = x + y_new * (new_mid_x + 1);

      fft_data_new[index_new] = fft_data[index_old] * factor;
    }
    if (width_in_ >= width_out_) {
      size_t index_old = width_in_ / 2 + y_old * (old_mid_x + 1);
      size_t index_new = width_out_ / 2 + y_new * (new_mid_x + 1);
      fft_data_new[index_new] = fft_data[index_old] * factor;
    }
  }

  fftwf_free(fft_data);

  fftwf_execute_dft_c2r(f_to_out_plan_,
                        reinterpret_cast<fftwf_complex*>(fft_data_new),
                        task.output);

  fftwf_free(fft_data_new);

  if (correct_window_ && window_function_ != WindowFunction::Rectangular &&
      !skip_window)
    UnapplyWindow(task.output);
}

void FFTResampler::MakeWindow(aocommon::UVector<float>& data,
                              size_t width) const {
  if (window_function_ == WindowFunction::Tukey)
    MakeTukeyWindow(data, width);
  else {
    data.resize(width);
    for (size_t x = 0; x != width; ++x)
      data[x] = WindowFunction::Evaluate(window_function_, width, x) + 1e-5;
  }
}

void FFTResampler::MakeTukeyWindow(aocommon::UVector<float>& data,
                                   size_t width) const {
  // Make a Tukey window, which consists of
  // left: a cosine going from 0 to 1
  // mid: all 1
  // right: a cosine going from 1 to 0
  data.resize(width);
  for (size_t x = 0; x != width; ++x) {
    // left part of Tukey window
    double x_sh = (0.5 + x) * 2;
    if (x_sh < width - tukey_inset_size_) {
      double pos = x_sh / (width - tukey_inset_size_);
      data[x] = (std::cos((pos + 1.0) * M_PI) + 1.0) * 0.5;
    } else if (x_sh < width + tukey_inset_size_) {
      data[x] = 1.0;
    } else {
      double pos =
          (x_sh - (width + tukey_inset_size_)) / (width - tukey_inset_size_);
      data[x] = (std::cos(pos * M_PI) + 1.0) * 0.5;
    }
  }
}

void FFTResampler::ApplyWindow(float* data) const {
  if (window_row_in_.empty()) {
    MakeWindow(window_row_in_, width_in_);
    MakeWindow(window_col_in_, height_in_);
    if (correct_window_) {
      aocommon::UVector<float> window_img_in(width_in_ * height_in_);
      float* ptr_in = window_img_in.data();
      for (size_t y = 0; y != height_in_; ++y) {
        for (size_t x = 0; x != width_in_; ++x) {
          *ptr_in = window_row_in_[x] * window_col_in_[y];
          ++ptr_in;
        }
      }

      window_out_.resize(width_out_ * height_out_);
      Task task;
      task.input = window_img_in.data();
      task.output = window_out_.data();
      RunSingle(task, true);
    }
  }
  for (size_t y = 0; y != height_in_; ++y) {
    for (size_t x = 0; x != width_in_; ++x) {
      *data *= window_row_in_[x] * window_col_in_[y];
      ++data;
    }
  }
}

void FFTResampler::UnapplyWindow(float* data) const {
  size_t n = width_out_ * height_out_;
  for (size_t i = 0; i != n; ++i) {
    data[i] /= window_out_[i];
  }
}
