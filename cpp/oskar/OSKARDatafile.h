#ifndef OSKAR_DATAFILE_H
#define OSKAR_DATAFILE_H

#include <string>
#include <memory>
#include <mutex>
#include <map>

#include <H5Cpp.h>

#include "OSKARDataset.h"

//! Oskar datafile interface
class Datafile {
 public:
  Datafile(const std::string& filename);

  std::shared_ptr<Dataset> get(const unsigned int freq);

 private:
  // Coeffs;
  std::map<unsigned int, std::shared_ptr<Dataset>> m_map;

  // HDF5
  std::string m_filename;
  std::unique_ptr<H5::H5File> m_h5_file;
  mutable std::mutex m_mutex;
};

#endif