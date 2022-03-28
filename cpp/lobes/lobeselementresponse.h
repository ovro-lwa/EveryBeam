// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_LOBES_ELEMENTRESPONSE_H
#define EVERYBEAM_LOBES_ELEMENTRESPONSE_H

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "../options.h"
#include "../fieldresponse.h"

#include <boost/optional.hpp>

#include <vector>
#include <memory>
#include <iostream>

namespace everybeam {

//! Implementation of the Lobes response model
class LOBESElementResponse : public FieldResponse {
 public:
  /**
   * @brief Construct a new LOBESElementResponse object
   *
   * @param name (LOFAR) station name, i.e. CS302LBA
   * @param options if options.coeff_path is non-empty it is used to find
   * coefficient files
   */
  LOBESElementResponse(const std::string& name, const Options& options);

  ElementResponseModel GetModel() const final override {
    return ElementResponseModel::kLOBES;
  }

  /**
   * @brief Stub override of the Response method, an element id
   * is needed to compute the element response.
   *
   * @param freq
   * @param theta
   * @param phi
   * @param response
   */
  aocommon::MC2x2 Response([[maybe_unused]] double freq,
                           [[maybe_unused]] double theta,
                           [[maybe_unused]] double phi) const final override {
    throw std::invalid_argument(
        "LOBESElementResponse::response needs an element_id");
  };

  /**
   * @brief Virtual implementation of Response method
   *
   * @param element_id ID of element
   * @param freq Frequency of the plane wave (Hz).
   * @param theta Angle wrt. z-axis (rad) (NOTE: parameter is just a stub to
   * match override)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad) (NOTE: parameter is
   * just a stub to match override)
   * @param result Pointer to 2x2 array of Jones matrix
   */
  aocommon::MC2x2 Response(int element_id, double freq, double theta,
                           double phi) const final override;

  /**
   * @brief Create LOBESElementResponse
   *
   * @param name Station name, e.g. CS302LBA
   * @return std::shared_ptr<LOBESElementResponse>
   */
  static std::shared_ptr<LOBESElementResponse> GetInstance(
      const std::string& name, const Options& options);

  /**
   * @brief Set field quantities (i.e. the basefunctions) for the LOBES element
   * response given the direction of interest, specified in (theta, phi) angles.
   * NOTE: this method overrides the "cached" basefunction_ member every time
   * when called.
   *
   * @param theta Angle wrt. z-axis (rad)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad)
   *
   * @warning When calling this function or @ref ClearFieldQuantities from
   * multiple threads the callers need to ensure the synchronisation.
   */
  virtual void SetFieldQuantities(double theta, double phi) final override {
    basefunctions_ = ComputeBaseFunctions(theta, phi);
  };

  /**
   * @brief Clear the cached basefunctions
   *
   * @warning When calling this function or @ref SetFieldQuantities from
   * multiple threads the callers need to ensure the synchronisation.
   */
  virtual void ClearFieldQuantities() final override {
    basefunctions_.reset();
  };

 private:
  // Typdef of BaseFunctions as Eigen::Array type
  typedef Eigen::Array<std::complex<double>, Eigen::Dynamic, 2> BaseFunctions;

  /** The cached version of the base functions. */
  mutable boost::optional<BaseFunctions> basefunctions_;

  // Find the closest frequency
  size_t FindFrequencyIdx(double f) const {
    auto is_closer = [f](int x, int y) { return abs(x - f) < abs(y - f); };
    auto result =
        std::min_element(frequencies_.begin(), frequencies_.end(), is_closer);
    return std::distance(frequencies_.begin(), result);
  }

  // Compute the base functions given theta and phi angles
  BaseFunctions ComputeBaseFunctions(double theta, double phi) const;

  // Store h5 coefficients in coefficients_
  Eigen::Tensor<std::complex<double>, 4, Eigen::RowMajor> coefficients_;

  // Store h5 frequencies in frequencies_
  std::vector<double> frequencies_;

  struct nms_t {
    int n, m, s;
  };
  std::vector<nms_t> nms_;
};
}  // namespace everybeam

#endif
