// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_OSKAR_DATASET_H_
#define EVERYBEAM_OSKAR_DATASET_H_

#include <complex>
#include <vector>

#include <H5Cpp.h>

#include <oskar_beam_utils.h>

namespace everybeam {
//! OSKAR dataset
class Dataset {
 public:
  /**
   * @brief Construct a new Dataset object given a h5 file and a
   * frequency
   *
   * @param h5_file H5 file (.h5)
   * @param freq Frequency to look for (Hz)
   */
  Dataset(H5::H5File& h5_file, const unsigned int freq);

  size_t GetNrElements() const { return nr_elements_; };
  size_t GetLMax() const { return l_max_; };

  const oskar::Double4C* GetAlphaPtr(const unsigned int element) const;

 private:
  // Using the OSKAR Double4C type ensures correct alignment.
  std::vector<oskar::Double4C> data_;
  size_t nr_elements_;
  size_t nr_coeffs_;
  size_t l_max_;
};
}  // namespace everybeam
#endif
