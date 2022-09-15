// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oskardatafile.h"

#include <iostream>
#include <tuple>

namespace everybeam {

Datafile::Datafile(const std::string& filename) {
  // Open file
  h5_file_.reset(new H5::H5File(filename, H5F_ACC_RDONLY));

  // Disable HDF5 error prints
  H5::Exception::dontPrint();
}

const Dataset& Datafile::Get(const unsigned int freq) {
  std::lock_guard<std::mutex> lock(mutex_);

  // Find dataset for frequency.
  auto entry = map_.find(freq);

  // If not found, read dataset.
  if (entry == map_.end()) {
    auto dataset = std::make_unique<Dataset>(*h5_file_, freq);
    std::tie(entry, std::ignore) = map_.insert({freq, std::move(dataset)});
  }

  // Return dataset.
  return *entry->second;
}
}  // namespace everybeam
