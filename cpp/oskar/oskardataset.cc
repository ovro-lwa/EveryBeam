// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <iostream>
#include <cassert>
#include <stdexcept>

#include <cstring>
#include <limits>
#include "oskardataset.h"

static constexpr unsigned kDatatsetRank = 3;

namespace everybeam {

Dataset::Dataset(H5::H5File& h5_file, const unsigned int freq) {
  H5::DataSet dataset;
  const int freq_mhz = (int)(freq / 1e6);
  try {
    // Try to find exact match
    const std::string dataset_name = std::to_string(freq_mhz);
    dataset = h5_file.openDataSet(dataset_name);
  } catch (H5::FileIException& e) {
    // Find nearest neighbor
    int dfreq = std::numeric_limits<int>::max();
    int freq_nearest = 0;
    for (size_t i = 0; i < h5_file.getNumObjs(); ++i) {
      const int h5_freq = stoi(h5_file.getObjnameByIdx(i));
      if (std::abs(h5_freq - freq_mhz) < dfreq) {
        dfreq = std::abs(h5_freq - freq_mhz);
        freq_nearest = h5_freq;
      }
    }
    std::cout << "Could not load dataset for frequency "
              << std::to_string(freq_mhz)
              << " MHz, using the nearest neighbor with frequency "
              << freq_nearest << " MHz instead" << std::endl;

    const std::string dataset_name = std::to_string(freq_nearest);
    try {
      dataset = h5_file.openDataSet(dataset_name);
    } catch (H5::FileIException& e) {
      throw std::runtime_error("Could not load dataset for frequency " +
                               dataset_name + "MHz.");
    }
  }

  // Read dataset dimensions
  H5::DataSpace dataspace = dataset.getSpace();
  const size_t rank = dataspace.getSimpleExtentNdims();
  assert(rank == kDatatsetRank);

  // Get dimensions
  std::vector<hsize_t> dims(rank);
  dataspace.getSimpleExtentDims(dims.data(), NULL);
  const size_t nr_elements = dims[0];
  const size_t nr_coeffs = dims[1];
  assert(dims[2] == 4);  // tetm*pol

// Coefficient data stored as:
// [nr_elements][nr_coefficients][4],
// with inner dimension:
// (x_te_re, x_te_im), (x_tm_re, x_tm_im),
// (y_te_re, y_te_im), (y_tm_re, y_tm_im)
#ifndef NDEBUG
  std::cout << "nr_elements: " << nr_elements << std::endl;
  std::cout << "nr_coeffs: " << nr_coeffs << std::endl;
#endif

  // Check total number of coefficients to find l_max
  const double l_max_d = sqrt(nr_coeffs + 1) - 1;
  const size_t l_max = (size_t)round(l_max_d);
#ifndef NDEBUG
  std::cout << "l_max: " << l_max << std::endl;
#endif

  // Sanity check
  assert(l_max * (l_max + 2) == nr_coeffs);

  // Set members
  nr_elements_ = nr_elements;
  nr_coeffs_ = nr_coeffs;
  l_max_ = l_max;

  // Read coefficients into data vector
  data_.resize(nr_elements * nr_coeffs * 4);
  assert(dims[0] * dims[1] * dims[2] == data_.size());
  H5::DataType data_type = dataset.getDataType();
  assert(data_type.getSize() == sizeof(std::complex<double>));
  dataset.read(data_.data(), data_type, dataspace);
}

size_t Dataset::GetIndex(const unsigned int element) const {
  return element * nr_coeffs_ * 4;
}

std::complex<double>* Dataset::GetAlphaPtr(const unsigned int element) {
  assert(element < GetNrElements());
  size_t index = GetIndex(element);
  return data_.data() + index;
}
}  // namespace everybeam
