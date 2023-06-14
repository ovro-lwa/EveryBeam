// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "atermresampler.h"

#include "../common/fftresampler.h"

#include <aocommon/fits/fitsreader.h>
#include <aocommon/imagecoordinates.h>

using aocommon::WindowFunction;

namespace everybeam {
namespace aterms {

ATermResampler::ATermResampler(
    const aocommon::CoordinateSystem& coordinate_system, size_t max_support)
    : coordinate_system_(coordinate_system),
      allocated_width_(max_support),
      allocated_height_(max_support),
      resampler_(),
      downsample_(true),
      window_(WindowFunction::RaisedHann),
      padding_(1.0),
      override_fits_phase_centre_(false),
      override_ra_(0.0),
      override_dec_(0.0) {}

ATermResampler::~ATermResampler() = default;

void ATermResampler::ReadAndResample(aocommon::FitsReader& reader,
                                     size_t file_index,
                                     aocommon::UVector<float>& scratch,
                                     aocommon::UVector<float>& output,
                                     double stretch_factor) {
  if (!resampler_) {
    resampler_.reset(new common::FFTResampler(
        allocated_width_, allocated_height_, coordinate_system_.width,
        coordinate_system_.height));
    if (window_ == WindowFunction::Tukey) {
      resampler_->SetTukeyWindow(double(allocated_width_) / padding_, false);
    } else {
      resampler_->SetWindowFunction(window_, true);
    }
  }

  if (downsample_) {
    reader.ReadIndex(output.data(), file_index);

    // First, the image is regridded on a smaller image that fits in the kernel
    // support allocated for the aterms
    regrid(reader, scratch.data(), output.data(), stretch_factor);

    // Now, the small image is enlarged so that it matches the kernel size
    resampler_->Resample(scratch.data(), output.data());
  } else {
    scratch.resize(reader.ImageWidth() * reader.ImageHeight());
    reader.ReadIndex(scratch.data(), file_index);

    regrid(reader, output.data(), scratch.data(), stretch_factor);
  }
}

size_t ATermResampler::ScratchBSize(const aocommon::FitsReader& reader) const {
  return std::max(coordinate_system_.width * coordinate_system_.height,
                  reader.ImageWidth() * reader.ImageHeight());
}

void ATermResampler::regrid(const aocommon::FitsReader& reader, float* dest,
                            const float* source, double stretch_factor) {
  using aocommon::ImageCoordinates;

  const size_t in_width = reader.ImageWidth();
  const size_t in_height = reader.ImageHeight();
  const double in_pixel_size_x = reader.PixelSizeX() / stretch_factor;
  const double in_pixel_size_y = reader.PixelSizeY() / stretch_factor;
  const double in_phase_centre_dl = reader.LShift();
  const double in_phase_centre_dm = reader.MShift();
  const double in_phase_centre_ra =
      override_fits_phase_centre_ ? override_ra_ : reader.PhaseCentreRA();
  const double in_phase_centre_dec =
      override_fits_phase_centre_ ? override_dec_ : reader.PhaseCentreDec();

  const size_t out_width =
      downsample_ ? allocated_width_ : coordinate_system_.width;
  const size_t out_height =
      downsample_ ? allocated_height_ : coordinate_system_.height;

  // The full size is regridded onto the 'Nyquist-sampled' image to remove
  // high-frequency components. aterm_dl/dm are the pixelsizes of the
  // Nyquist-sampled image.
  const double aterm_dl =
      coordinate_system_.dl * coordinate_system_.width / out_width;
  const double aterm_dm =
      coordinate_system_.dm * coordinate_system_.height / out_height;
  /**
   * If phase centra of input and output are the same, i.e. they have the same
   * tangential plane, a few calculations can be saved.
   */
  const bool same_plane = in_phase_centre_ra == coordinate_system_.ra &&
                          in_phase_centre_dec == coordinate_system_.dec;

  size_t index = 0;
  for (size_t y = 0; y != out_width; ++y) {
    for (size_t x = 0; x != out_width; ++x) {
      double l, m;
      ImageCoordinates::XYToLM(x, y, aterm_dl, aterm_dm, out_width, out_width,
                               l, m);
      l += coordinate_system_.l_shift;
      m += coordinate_system_.m_shift;
      if (!same_plane) {
        double pixra, pixdec;
        ImageCoordinates::LMToRaDec(l, m, coordinate_system_.ra,
                                    coordinate_system_.dec, pixra, pixdec);
        ImageCoordinates::RaDecToLM(pixra, pixdec, in_phase_centre_ra,
                                    in_phase_centre_dec, l, m);
      }
      l -= in_phase_centre_dl;
      m -= in_phase_centre_dm;
      int inX, inY;
      ImageCoordinates::LMToXY(l, m, in_pixel_size_x, in_pixel_size_y, in_width,
                               in_height, inX, inY);
      if (inX < 0 || inY < 0 || inX >= int(in_width) || inY >= int(in_height)) {
        dest[index] = 0;
      } else {
        dest[index] = source[inX + inY * in_width];
      }
      ++index;
    }
  }
}

}  // namespace aterms
}  // namespace everybeam
