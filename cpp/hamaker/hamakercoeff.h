// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_HAMAKER_COEFF_H_
#define EVERYBEAM_HAMAKER_COEFF_H_

#include <iostream>
#include <string>
#include <complex>
#include <cassert>
#include <vector>
#include <cstring>

#include <H5Cpp.h>

namespace everybeam {
//! Hamaker coefficients
class HamakerCoefficients {
 public:
  //! Default constructor
  HamakerCoefficients();

  //! Constructor for reading coeff from file
  HamakerCoefficients(const std::string& filename);

  //! Constructor for writing coeff to file
  HamakerCoefficients(const double freq_center, const double freq_range,
                      const unsigned int nHarmonics,
                      const unsigned int nPowerTheta,
                      const unsigned int nPowerFreq);

  /**
   * @brief Set Hamaker coefficients
   *
   * @param n
   * @param t
   * @param f
   * @param value
   */
  void SetCoefficients(
      unsigned int n, unsigned int t, unsigned int f,
      std::pair<std::complex<double>, std::complex<double>> value);

  void SetCoefficients(const std::complex<double>* coeff);

  void SetCoefficients(const std::vector<std::complex<double>> coeff);

  // Get
  size_t GetNrCoefficients() const;

  double GetFreqCenter() const { return freq_center_; }

  double GetFreqRange() const { return freq_range_; }

  unsigned int Get_nHarmonics() const { return nHarmonics_; }

  unsigned int Get_nPowerTheta() const { return nPowerTheta_; }

  unsigned int Get_nPowerFreq() const { return nPowerFreq_; }

  void GetCoefficient(
      unsigned int n, unsigned int t, unsigned int f,
      std::pair<std::complex<double>, std::complex<double>>& value) const;

  // HDF5 I/O
  void ReadCoefficients(const std::string& filename);

  void WriteCoefficients(const std::string& filename);

  // Debugging
  void PrintCoefficients();

 private:
  // Methods
  size_t GetIndex(unsigned int n, unsigned int t, unsigned int f) const;

  // Parameters
  double freq_center_;
  double freq_range_;
  unsigned int nHarmonics_;
  unsigned int nPowerTheta_;
  unsigned int nPowerFreq_;
  const unsigned int nInner_ = 2;

  // Data
  std::vector<std::complex<double>> coeff_;

  // HDF5
  std::string dataset_name_ = "coeff";
  const unsigned int dataset_rank_ = 4;
};
}  // namespace everybeam
#endif
