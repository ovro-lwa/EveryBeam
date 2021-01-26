// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_H5ATERM_H
#define EVERYBEAM_ATERMS_H5ATERM_H

#include "atermbase.h"
#include "cache.h"
#include "../coords/coordutils.h"

#include <complex>
#include <map>
#include <memory>
#include <vector>

#include <aocommon/uvector.h>
#include <schaapcommon/h5parm/h5parm.h>
#include <cassert>
#include <numeric>

namespace everybeam {
namespace aterms {

/**
 * @brief Convenience class for efficiently evaluating
 * a Lagrange binomial (i.e. polynomial in 2D), using
 * Horner's method.
 *
 */
class LagrangePolynomial {
 public:
  /**
   * @brief Construct a new Lagrange Polynomial object
   *
   * @param nr_coeffs
   */
  LagrangePolynomial(size_t nr_coeffs) : nr_coeffs_(nr_coeffs) {
    order_ = ComputeOrder(nr_coeffs_);
  }

  /**
   * @brief Evaluate binomial on x,y coordinate
   *
   * @param x x-coordinate
   * @param y y-coordinate
   * @param coeffs Coefficients of polynomial. Ordering as in Pascal's triangle
   * is assumed
   * @return float
   */
  float Evaluate(double x, double y, const std::vector<float>& coeffs) {
    assert(coeffs.size() == nr_coeffs_);
    std::vector<float> y_coeffs(order_ + 1);
    for (int i = static_cast<int>(order_); i >= 0; --i) {
      std::vector<size_t> coeff_indices = GetRDiagonalIndices(i);
      std::vector<float> coeffs_tmp;
      for (auto idx : coeff_indices) {
        coeffs_tmp.push_back(coeffs[idx]);
      }
      // Flip ordering of coefficients
      std::reverse(coeffs_tmp.begin(), coeffs_tmp.end());
      y_coeffs[i] = HornersMethod(y, coeffs_tmp);
    }
    // Flip ordering
    std::reverse(y_coeffs.begin(), y_coeffs.end());
    float result = HornersMethod(x, y_coeffs);
    return result;
  }

  // Get indices of right diagonal in Pascal's triangle,
  // given the number of the diagonal
  std::vector<size_t> GetRDiagonalIndices(size_t diag_nr) {
    if (diag_nr > order_) {
      throw std::runtime_error(
          "Diagonal number should be smaller/equal to polynomial order");
    }
    // Offset per diaonal equals nr of terms that correspond to nr of terms in
    // polynomial of one order lower
    std::vector<size_t> indices = {
        (diag_nr == 0) ? diag_nr : ComputeNrCoeffs(diag_nr - 1)};
    for (size_t i = diag_nr + 1; i <= order_; ++i) {
      indices.push_back(indices.back() + i + 1);
    }
    return indices;
  }

  /**
   * @brief Compute polynomial order, given the total number
   * of coefficients
   *
   * @param nr_coeffs Number of coefficients
   * @return size_t Polynomial order
   */
  static size_t ComputeOrder(size_t nr_coeffs) {
    size_t degree = 0;
    size_t term_counter = 1;
    while (term_counter < nr_coeffs) {
      degree++;
      term_counter += (degree + 1);
    }

    // Computed number of coeffs and input value should match exactly
    if (term_counter != nr_coeffs) {
      throw std::runtime_error(
          "Computed number of coefficients and input number of coefficients do "
          "not match exactly. Probably the input value is inconsistent.");
    }
    return degree;
  }

  /**
   * @brief Compute number of coeffs, given the polynomial order
   *
   * @param order polynomial order
   * @return size_t number of terms
   */
  static size_t ComputeNrCoeffs(size_t order) {
    size_t nr_coeffs = 0;
    for (size_t i = 0; i <= order; ++i) {
      nr_coeffs += (i + 1);
    }
    return nr_coeffs;
  }

  size_t GetOrder() const { return order_; }
  size_t GetNrCoeffs() const { return nr_coeffs_; }

 private:
  /**
   * @brief Evaluate (1D) polynomial, using Horner's method
   *
   * @param coord Evaluation coordinate
   * @param coeffs Coefficients, need to be ordered as in c_0 + c_1 x + ... +
   * c_n x^n
   * @return float Result
   */
  static float HornersMethod(double coord, const std::vector<float>& coeffs) {
    float result = 0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
      result = coeffs[i] + result * coord;
    }
    return result;
  }

  size_t nr_coeffs_, order_;
};

/**
 * Class that reads in H5Parm coefficient files and evaluates the
 * underlying polynomial on a prescribed image.
 * The H5Parm file(s) are supposed to have an "amplitude_coefficients"
 * and a "phase_coefficients" solution table, where each solution table
 * has at least the following axes ("ant", "time", "dir"). The polynomial
 * coefficients are stored along the "dir" axis
 */
class H5ParmATerm : public ATermBase {
 public:
  H5ParmATerm(const std::vector<std::string>& station_names_ms,
              const coords::CoordinateSystem& coordinate_system);

  ~H5ParmATerm(){};

  /**
   * @brief Read h5parm files given a vector of paths
   *
   * @param filenames
   */
  void Open(const std::vector<std::string>& filenames);

  /**
   * @brief
   *
   * @param buffer Buffer
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s)
   * @param frequency Freq (Hz) - not used at the moment
   * @param field_id Irrelevant for h5parm aterms
   * @param uvw_in_m Irrelevant for h5parm aterms
   * @return true Results are updated
   * @return false No need to update the result, cached result can be used
   */
  bool Calculate(std::complex<float>* buffer, double time, double frequency,
                 size_t field_id, const double* uvw_in_m) override final;

  /**
   * @brief Set the update interval
   *
   * @param update_interval Update interval (in s)
   */
  void SetUpdateInterval(double update_interval) {
    update_interval_ = update_interval;
  }

  // Get average update time, fixed value for h5parm aterm
  double AverageUpdateTime() const override final { return update_interval_; }

 private:
  // Expand complex exponential from amplitude and phase
  std::complex<float> ExpandComplexExp(const std::string& station_name,
                                       hsize_t ampl_tindex,
                                       hsize_t phase_tindex, double freq,
                                       double l, double m,
                                       bool recalculate_ampl,
                                       bool recalculate_phase, size_t offset);

  // Read coefficients from solution tab, for given
  // time index(frequency not relevant, as yet)
  static void ReadCoeffs(schaapcommon::h5parm::SolTab& soltab,
                         const std::string& station_name,
                         std::vector<float>& coeffs, hsize_t tindex, double);

  std::vector<schaapcommon::h5parm::SolTab> amplitude_soltab_;
  std::vector<schaapcommon::h5parm::SolTab> phase_soltab_;
  const std::vector<std::string> station_names_ms_;

  // Store polynomial information
  std::unique_ptr<LagrangePolynomial> ampl_polynomial_;
  std::unique_ptr<LagrangePolynomial> phase_polynomial_;
  coords::CoordinateSystem coordinate_system_;

  // Top level (i.e. ATermConfig) caching
  double update_interval_;
  double last_aterm_update_;

  // Amplitude and phase caching
  hsize_t last_ampl_index_;
  hsize_t last_phase_index_;
  aocommon::UVector<float> amplitude_cache_;
  aocommon::UVector<float> phase_cache_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
