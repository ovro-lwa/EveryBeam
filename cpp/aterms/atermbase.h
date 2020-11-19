// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_ATERM_BASE_H
#define EVERYBEAM_ATERMS_ATERM_BASE_H

#include <complex>
#include <string>

namespace everybeam {
namespace aterms {

class ATermBase {
 public:
  ATermBase() : save_aterms_(false) {}

  virtual ~ATermBase() {}

  /**
   * Calculate the a-terms for the given time and frequency, for all stations.
   * @param buffer A buffer of size 4 x width x height x nstation, in that
   * order.
   * @param time The time corresponding to the currently gridded visibilities to
   * which the aterm will be applied.
   * @param frequency Frequency of currently gridded visibilities.
   * @param fieldId Field that these visibilities belong to. Different fields
   * might requires different beams.
   * @param uvwInM The UVW of the antennas, referenced to antenna 0.
   * @returns @c True when new aterms are calculated. If these aterms are the
   * same as for the previous call to Calculate(), @c false can be returned and
   * the output buffer does not need to be updated. The gridder will then make
   * sure to use the previous aterms, and not reserve extra memory for it etc.
   */
  virtual bool Calculate(std::complex<float>* buffer, double time,
                         double frequency, size_t fieldId,
                         const double* uvwInM) = 0;

  virtual double AverageUpdateTime() const = 0;

  void StoreATermsEigenvalues(const std::string& filename,
                              const std::complex<float>* buffer,
                              size_t n_stations, size_t width, size_t height);

  void StoreATermsReal(const std::string& filename,
                       const std::complex<float>* buffer, size_t n_stations,
                       size_t width, size_t height);

  /**
   * Set whether a fits image with the a-terms should be written to disk
   * every time they are calculated.
   * @param save_aterms Fits images are saved when set to true.
   */
  void SetSaveATerms(bool save_aterms, const std::string& prefix) {
    save_aterms_ = save_aterms;
    prefix_ = prefix;
  }

 protected:
  void SaveATermsIfNecessary(const std::complex<float>* buffer,
                             size_t n_stations, size_t width, size_t height) {
    if (save_aterms_) {
      static int index = 0;
      std::ostringstream f;
      StoreATermsEigenvalues(
          prefix_ + "-aterm-ev" + std::to_string(index) + ".fits", buffer,
          n_stations, width, height);
      StoreATermsReal(
          prefix_ + "-aterm-realxx" + std::to_string(index) + ".fits", buffer,
          n_stations, width, height);
      ++index;
    }
  }

 private:
  bool save_aterms_;
  std::string prefix_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
