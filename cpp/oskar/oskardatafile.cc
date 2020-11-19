// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oskardatafile.h"

#include <iostream>

Datafile::Datafile(const std::string& filename) {
  // Open file
  std::cout << "read oskar datafile: " << filename << std::endl;
  h5_file_.reset(new H5::H5File(filename, H5F_ACC_RDONLY));

  // Disable HDF5 error prints
  H5::Exception::dontPrint();
};

std::shared_ptr<Dataset> Datafile::Get(const unsigned int freq) {
  std::lock_guard<std::mutex> lock(mutex_);

  // Find dataset for frequency
  auto entry = map_.find(freq);

  // If found, retrieve pointer to dataset
  if (entry != map_.end()) {
    return entry->second;
  }

  // Read and return dataset
  std::shared_ptr<Dataset> dataset_ptr;
  dataset_ptr.reset(new Dataset(*h5_file_, freq));
  map_.insert({freq, dataset_ptr});
  return dataset_ptr;
}
