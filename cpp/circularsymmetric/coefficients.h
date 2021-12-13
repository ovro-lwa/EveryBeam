// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_CIRCULARSYMMETRIC_COEFFICIENTS_H_
#define EVERYBEAM_CIRCULARSYMMETRIC_COEFFICIENTS_H_

#include <array>
#include <map>
#include <string>

#include <aocommon/uvector.h>

namespace everybeam {
namespace circularsymmetric {

/**
 * Holds information about the polynomial coefficients of the directional
 * response of a telescope. It also contains some further information to
 * construct the response. The coefficients make up a polynomial that construct
 * a circularly symmetric beam. Only even coefficients are specified, i.e. the
 * function to evaluate is:
 *
 * response(d) = a_0 + a_1 d^2 + a_2 d^4 + ...
 *
 * The class @ref VoltagePattern can calculate the response from a coefficient
 * class.
 */
class Coefficients {
 public:
  virtual ~Coefficients() = default;
  /**
   * List of @c n_coef x @c n_freq polynomial coefficients. If this list
   * contains multiple coefficients for multiple frequencies, the @ref
   * VoltagePattern class will interpolate them linearly. If interpolation is
   * not desired/necessary, the frequency parameter is used to select the right
   * coefficients and return a list of size @c n_coef. In the latter case,
   * @ref GetFrequencies() returns only one value.
   * @param frequency Frequency in Hz.
   */
  virtual aocommon::UVector<double> GetCoefficients(double frequency) const = 0;
  /**
   * List of frequencies in units of Hz, corresponding to the coefficients
   * returned by
   * @ref GetCoefficients(). The frequencies are used for interpolation. If
   * interpolation is not necessary / desirable, this function returns only one
   * element.
   * @param frequency Frequency in Hz.
   */
  virtual aocommon::UVector<double> GetFrequencies(double frequency) const = 0;
  /**
   * @returns How far out the beam is evaluated in the one dimensional (radial)
   * direction, in units of arc minutes. This does not influence the shape of
   * the beam, only the furthest distance able to calculate.
   */
  virtual double MaxRadiusInArcMin() const = 0;
  /**
   * @returns A reference frequency in units of Hz. Most beams are expressed in
   * arcmin per GHz and should therefore return 1e9, but some beams are
   * calculated with a different reference. The beam is stretched or shrunk when
   * the reference frequency is decreased or increased, respectively.
   */
  virtual double ReferenceFrequency() const = 0;
  /**
   * True if the coefficients make up an inverted beam. An inverted beam is
   * one that converts apparent values to intrinsic values. In the inverted
   * case, the second coefficient is positive instead of negative.
   */
  virtual bool AreInverted() const = 0;
};
}  // namespace circularsymmetric
}  // namespace everybeam

#endif  // EVERYBEAM_CIRCULARSYMMETRIC_COEFFICIENTS_H_
