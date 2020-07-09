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
  size_t get_nr_elements() const { return m_nr_elements; };
  size_t GetLMax() const { return m_l_max; };

  std::complex<double>* GetAlphaPtr(const unsigned int element);

 private:
  // Methods
  size_t GetIndex(const unsigned int element) const;

  // Constants
  const unsigned int dataset_rank_ = 3;

  // Members
  std::vector<std::complex<double>> m_data;
  unsigned int m_nr_elements;
  unsigned int m_nr_coeffs;
  unsigned int m_l_max;
};

#endif