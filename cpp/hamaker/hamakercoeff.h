#ifndef HAMAKER_COEFF_H
#define HAMAKER_COEFF_H

#include <iostream>
#include <string>
#include <complex>
#include <cassert>
#include <vector>
#include <cstring>

#include <H5Cpp.h>

//! Hamaker coefficients
class HamakerCoefficients {
 public:
  //! Default constructor
  HamakerCoefficients();

  //! Constructor for reading coeff from file
  HamakerCoefficients(std::string& filename);

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
      const unsigned int n, const unsigned int t, const unsigned int f,
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

  std::pair<std::complex<double>, std::complex<double>> GetCoefficient(
      const unsigned int n, const unsigned int t, const unsigned int f);

  // HDF5 I/O
  void ReadCoefficients(std::string& filename);

  void WriteCoefficients(std::string& filename);

  // Debugging
  void PrintCoefficients();

 private:
  // Methods
  size_t GetIndex(const unsigned int n, const unsigned int t,
                  const unsigned int f);

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

#endif
