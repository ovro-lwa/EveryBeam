#ifndef EVERYBEAM_ELEMENT_HAMAKER_H
#define EVERYBEAM_ELEMENT_HAMAKER_H

#include <complex>
#include <memory>

#include "antenna.h"
#include "element.h"
#include "elementresponse.h"
#include "common/types.h"

namespace everybeam {

/**
 * @brief Elementary antenna, optimized for LOFAR Hamaker model. Derived from
 * the Element class.
 *
 */
class ElementHamaker : public Element {
 public:
  /**
   * @brief Construct a new Element object
   *
   * @param coordinate_system (antenna) CoordinateSystem
   * @param element_response ElementResponseModel
   * @param id
   */
  ElementHamaker(const CoordinateSystem &coordinate_system,
                 ElementResponse::Ptr element_response, int id)
      : Element(coordinate_system, element_response, id) {}

  Antenna::Ptr Clone() const override;

  virtual matrix22c_t Response(real_t time, real_t freq,
                               const vector3r_t &direction,
                               const Options &options) override {
    // The only transform that is needed is hard-coded in LocalResponse
    matrix22c_t response = LocalResponse(time, freq, direction, options);
    return response;
  }

  virtual matrix22c_t LocalResponse(real_t time, real_t freq,
                                    const vector3r_t &direction, size_t id,
                                    const Options &options) const override;

 private:
  virtual matrix22c_t LocalResponse(
      real_t time, real_t freq, const vector3r_t &direction,
      const Options &options) const final override {
    return LocalResponse(time, freq, direction, id_, options);
  };
};
}  // namespace everybeam

#endif
