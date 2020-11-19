// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "atermbase.h"

#include <aocommon/fits/fitswriter.h>
#include <aocommon/matrix2x2.h>
#include <aocommon/uvector.h>

using aocommon::FitsWriter;
using aocommon::Matrix2x2;
using everybeam::aterms::ATermBase;

void ATermBase::StoreATermsEigenvalues(const std::string& filename,
                                       const std::complex<float>* buffer,
                                       size_t n_stations, size_t width,
                                       size_t height) {
  const size_t ny = std::floor(std::sqrt(n_stations));
  const size_t nx = (n_stations + ny - 1) / ny;
  aocommon::UVector<double> img(nx * ny * width * height, 0.0);
  for (size_t ant = 0; ant != n_stations; ++ant) {
    size_t x_corner = (ant % nx) * width, y_corner = (ant / nx) * height;
    for (size_t y = 0; y != height; ++y) {
      for (size_t x = 0; x != width; ++x) {
        std::complex<float> e1, e2;
        Matrix2x2::EigenValues(buffer + (width * (ant * height + y) + x) * 4,
                               e1, e2);
        double val = std::max(std::abs(e1), std::abs(e2));
        img[(y_corner + y) * width * nx + x + x_corner] = val;
      }
    }
  }
  FitsWriter writer;
  writer.SetImageDimensions(nx * width, ny * height);
  writer.Write(filename, img.data());
}

void ATermBase::StoreATermsReal(const std::string& filename,
                                const std::complex<float>* buffer,
                                size_t n_stations, size_t width,
                                size_t height) {
  size_t ny = std::floor(std::sqrt(n_stations)),
         nx = (n_stations + ny - 1) / ny;
  aocommon::UVector<double> img(nx * ny * width * height, 0.0);
  for (size_t ant = 0; ant != n_stations; ++ant) {
    size_t x_corner = (ant % nx) * width, y_corner = (ant / nx) * height;
    for (size_t y = 0; y != height; ++y) {
      for (size_t x = 0; x != width; ++x) {
        std::complex<float> xx =
            (buffer + (width * (ant * height + y) + x) * 4)[0];
        img[(y_corner + y) * width * nx + x + x_corner] = xx.real();
      }
    }
  }
  FitsWriter writer;
  writer.SetImageDimensions(nx * width, ny * height);
  writer.Write(filename, img.data());
}
