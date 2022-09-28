// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ELEMENTRESPONSE_H
#define EVERYBEAM_ELEMENTRESPONSE_H

#include <complex>
#include <filesystem>
#include <memory>
#include <ostream>
#include <vector>

#include <aocommon/matrix2x2.h>

#include "common/types.h"

namespace everybeam {

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
 * This class uses std::enable_shared_from_this since
 * LobesElementResponseFixedDirection stores a pointer to the wrapped class,
 * which is already owned elsewhere (typically by Station).
 */
class ElementResponse : public std::enable_shared_from_this<ElementResponse> {
 public:
  virtual ~ElementResponse() {}

  virtual ElementResponseModel GetModel() const = 0;

  /**
   * Create an element response object for a fixed direction.
   * This function allows reusing direction-specific values in
   * the newly created ElementResponse object.
   *
   * This default implementation creates an ElementResponseFixedDirection,
   * which fixates the theta and phi values for all Response() calls.
   * @param direction Direction of arrival (ITRF, m).
   * @return A new ElementResponse object, which should be used instead of the
   *         current ElementResponse object.
   */
  [[nodiscard]] virtual std::shared_ptr<ElementResponse> FixateDirection(
      const vector3r_t& direction) const;

  /**
   * @brief Virtual implementation of Response method
   *
   * @param cache Cached data from CacheDirection().
   * @param freq Frequency of the plane wave (Hz).
   * @param theta Angle wrt. z-axis (rad)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad)
   * @return A 2x2 array containing the Jones matrix.
   */
  virtual aocommon::MC2x2 Response(double freq, double theta,
                                   double phi) const = 0;

  /**
   * @brief Virtual implementation of Response method
   *
   * @param cache Cached data from CacheDirection().
   * @param element_id ID of element
   * @param freq Frequency of the plane wave (Hz).
   * @param theta Angle wrt. z-axis (rad)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad)
   * @return A 2x2 array containing the Jones matrix.
   */
  virtual aocommon::MC2x2 Response([[maybe_unused]] int element_id, double freq,
                                   double theta, double phi) const {
    return Response(freq, theta, phi);
  }

  static std::shared_ptr<const ElementResponse> GetInstance(
      ElementResponseModel model, const std::string& name,
      const Options& options);

 protected:
  /**
   * Get the path to an EveryBeam data file or directory.
   * @param relative_path A path relative to the EveryBeam data directory.
   * @return The full path of the data file or data directory.
   */
  static std::filesystem::path GetPath(
      const std::filesystem::path& relative_path);
};
}  // namespace everybeam
#endif
