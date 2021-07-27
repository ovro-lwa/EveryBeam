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

  std::shared_ptr<Dataset> Get(const unsigned int freq);

 private:
  // Coeffs;
  std::map<unsigned int, std::shared_ptr<Dataset>> map_;

  // HDF5
  std::string filename_;
  std::unique_ptr<H5::H5File> h5_file_;
  mutable std::mutex mutex_;
};
}  // namespace everybeam

#endif