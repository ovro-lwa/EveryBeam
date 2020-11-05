#ifndef EVERYBEAM_FIELDRESPONSE_H
#define EVERYBEAM_FIELDRESPONSE_H

#include "elementresponse.h"

namespace everybeam {

/**
 * @brief Virtual class that inherits from the ElementResponse class but can be
 * used to precompute quantities at "field" level. Currently, the
 * LOBESElementResponse class inherits from this class
 *
 */
class FieldResponse : public ElementResponse {
 public:
  /**
   * @brief Virtual method that can be used to set field quantities based
   * on the direction of interest, specified in (theta, phi) angles
   *
   * @param theta Angle wrt. z-axis (rad)
   * @param phi Angle in the xy-plane wrt. x-axis  (rad)
   */
  virtual void SetFieldQuantities(double theta, double phi) = 0;

  /**
   * @brief Clear the cached field quantities
   *
   */
  virtual void ClearFieldQuantities() = 0;
};
}  // namespace everybeam
#endif
