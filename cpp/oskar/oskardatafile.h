// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_OSKAR_DATAFILE_H_
#define EVERYBEAM_OSKAR_DATAFILE_H_

#include <string>
#include <memory>
#include <mutex>
#include <map>

#include <H5Cpp.h>

#include "oskardataset.h"

namespace everybeam {

//! Oskar datafile interface
class Datafile {
 public:
  Datafile(const std::string& filename);

  /**
   * Get the Dataset for a frequency.
   *
   * The first request for a Dataset for a frequency loads the data from
   * from the HDF5 file. Successive requests reuse the previously read Dataset.
   *
   * @param freq Frequency identifier.
   * @return A reference to the Dataset. It remains valid as long as the
   * Datafile that returned the Dateset is valid.
   */
  const Dataset& Get(const unsigned int freq);

 private:
  // Maps frequency identifiers to Datasets.
  std::map<unsigned int, std::unique_ptr<Dataset>> map_;

  // HDF5
  std::string filename_;
  std::unique_ptr<H5::H5File> h5_file_;
  // Protects concurrent access to map_.
  std::mutex mutex_;
};
}  // namespace everybeam

#endif