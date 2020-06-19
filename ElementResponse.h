#ifndef EVERYBEAM_ELEMENTRESPONSE_H
#define EVERYBEAM_ELEMENTRESPONSE_H

#include <complex>
#include <ostream>

#include "MutablePtr.h"

namespace everybeam {

enum ElementResponseModel {
  Unknown,
  Hamaker,
  LOBES,
  OSKARDipole,
  OSKARSphericalWave
};

std::ostream& operator<<(std::ostream& os, ElementResponseModel model);

/**
 * @brief Virtual class for the element response model. All the
 * (antenna/element) response models inherit from this class.
 *
 */
class ElementResponse {
 public:
  typedef MutablePtr<ElementResponse>
      Ptr;  //!< Pointer to ElementResponse object

  /**
   * @brief Virtual implementation of response method
   *
   * @param freq Frequency of the plane wave (Hz).
   * @param theta Angle wrt. z-axis (rad)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad)
   * @param result Pointer to 2x2 array of Jones matrix
   */
  virtual void response(double freq, double theta, double phi,
                        std::complex<double> (&result)[2][2]) const = 0;

  /**
   * @brief Virtual implementation of response method
   *
   * @param element_id ID of element
   * @param freq Frequency of the plane wave (Hz).
   * @param theta Angle wrt. z-axis (rad)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad)
   * @param result Pointer to 2x2 array of Jones matrix
   */
  virtual void response(int element_id, double freq, double theta, double phi,
                        std::complex<double> (&result)[2][2]) const {
    response(freq, theta, phi, result);
  }
};
}  // namespace everybeam
#endif
