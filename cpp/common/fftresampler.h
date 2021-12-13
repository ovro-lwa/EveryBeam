// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_COMMON_FFT_RESAMPLE_H_
#define EVERYBEAM_COMMON_FFT_RESAMPLE_H_

#include <aocommon/windowfunction.h>
#include <aocommon/uvector.h>

#include <vector>
#include <thread>

#include <fftw3.h>

namespace everybeam {
namespace common {

/**
 * @brief FFT resampling from (coarse) input grid to (high resolution) output
 * grid
 *
 */
class FFTResampler {
 private:
  struct Task {
    float *input, *output;
  };

 public:
  /**
   * @brief Construct a new FFTResampler object
   *
   * @param width_in input image width (int)
   * @param height_in input image height (int)
   * @param width_out output image width (int)
   * @param height_out output image height (int)
   */
  FFTResampler(size_t width_in, size_t height_in, size_t width_out,
               size_t height_out);

  ~FFTResampler();

  /**
   * @brief Do the FFT resampling
   *
   * @param input Input image buffer
   * @param output Output image buffer
   */
  void Resample(float* input, float* output) {
    if (width_in_ == width_out_ && height_in_ == height_out_) {
      std::copy_n(input, width_in_ * height_in_, output);
    } else {
      Task task;
      task.input = input;
      task.output = output;
      RunSingle(task, false);
    }
  }

  /**
   * Only to be used with SingleFT (it makes resampling thread unsafe!)
   */
  void SetTukeyWindow(double inset_size, bool correct_window) {
    window_function_ = aocommon::WindowFunction::Tukey;
    tukey_inset_size_ = inset_size;
    correct_window_ = correct_window;
    window_row_in_.clear();
    window_col_in_.clear();
    window_out_.clear();
  }

  void SetWindowFunction(aocommon::WindowFunction::Type window,
                         bool correct_window) {
    window_function_ = window;
    correct_window_ = correct_window;
    window_row_in_.clear();
    window_col_in_.clear();
    window_out_.clear();
  }

 private:
  void RunSingle(const Task& task, bool skip_window) const;
  void ApplyWindow(float* data) const;
  void UnapplyWindow(float* data) const;
  void MakeWindow(aocommon::UVector<float>& data, size_t width) const;
  void MakeTukeyWindow(aocommon::UVector<float>& data, size_t width) const;

  size_t width_in_, height_in_;
  size_t width_out_, height_out_;
  size_t fft_width_, fft_height_;
  aocommon::WindowFunction::Type window_function_;
  double tukey_inset_size_;
  mutable aocommon::UVector<float> window_row_in_;
  mutable aocommon::UVector<float> window_col_in_;
  mutable aocommon::UVector<float> window_out_;
  bool correct_window_;

  fftwf_plan in_to_f_plan_, f_to_out_plan_;
};
}  // namespace common
}  // namespace everybeam
#endif  // EVERYBEAM_COMMON_FFT_RESAMPLE_H_
