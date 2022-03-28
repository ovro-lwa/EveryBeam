// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ELEMENTRESPONSE_H
#define EVERYBEAM_ELEMENTRESPONSE_H

#include <complex>
#include <ostream>
#include <aocommon/matrix2x2.h>

#include "common/mutable_ptr.h"

namespace everybeam {

namespace common {
template <typename T>
class MutablePtr;
}

struct Options;

enum class ElementResponseModel {
  /// The default will select the default element response model
  /// based on the telescope: e.g. LOFAR will use kHamaker,
  /// OSKAR will select kOSKARSphericalWave.
  ///
  /// @note The values are part of the public ABI.
  //
  // kHamakerLba is specifically introduced for AARTFAAC observations, in which
  // case the type HBA/LBA cannot be inferred from the station name
  kDefault = 0,
  kHamaker = 1,
  kHamakerLba = 2,
  kLOBES = 3,
  kOSKARDipole = 4,
  kOSKARSphericalWave = 5,
  kSkaMidAnalytical = 6,
  kAartfaacInner = 7,
  kAartfaacOuter = 8
};

std::ostream& operator<<(std::ostream& os, ElementResponseModel model);

ElementResponseModel ElementResponseModelFromString(
    const std::string& element_response);

/**
 * @brief Abstract class for the element response model. All the
 * (antenna/element) response models inherit from this class.
 *
 */
class ElementResponse {
 public:
  virtual ~ElementResponse() {}

  typedef common::MutablePtr<ElementResponse>
      Ptr;  //!< Pointer to ElementResponse object

  virtual ElementResponseModel GetModel() const = 0;

  /**
   * @brief Virtual implementation of Response method
   *
   * @param freq Frequency of the plane wave (Hz).
   * @param theta Angle wrt. z-axis (rad)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad)
   * @param result Pointer to 2x2 array of Jones matrix
   */
  virtual aocommon::MC2x2 Response(double freq, double theta,
                                   double phi) const = 0;

  /**
   * @brief Virtual implementation of Response method
   *
   * @param element_id ID of element
   * @param freq Frequency of the plane wave (Hz).
   * @param theta Angle wrt. z-axis (rad)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad)
   * @param result Pointer to 2x2 array of Jones matrix
   */
  virtual aocommon::MC2x2 Response([[maybe_unused]] int element_id, double freq,
                                   double theta, double phi) const {
    return Response(freq, theta, phi);
  }

  static std::shared_ptr<ElementResponse> GetInstance(
      ElementResponseModel model, const std::string& name,
      const Options& options);
};
}  // namespace everybeam
#endif
