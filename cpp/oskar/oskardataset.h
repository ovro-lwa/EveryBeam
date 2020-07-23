#ifndef OSKAR_DATASET_H
#define OSKAR_DATASET_H

#include <complex>
#include <vector>

#include <H5Cpp.h>

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

  // Get
  size_t GetNrElements() const { return nr_elements_; };
  size_t GetLMax() const { return l_max_; };

  std::complex<double>* GetAlphaPtr(const unsigned int element);

 private:
  // Methods
  size_t GetIndex(const unsigned int element) const;

  // Constants
  const unsigned int dataset_rank_ = 3;

  // Members
  std::vector<std::complex<double>> data_;
  unsigned int nr_elements_;
  unsigned int nr_coeffs_;
  unsigned int l_max_;
};

#endif