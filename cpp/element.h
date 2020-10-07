#ifndef EVERYBEAM_ELEMENT_H
#define EVERYBEAM_ELEMENT_H

#include <complex>
#include <memory>

#include "antenna.h"
#include "elementresponse.h"
#include "common/types.h"

namespace everybeam {

/**
 * @brief Elementary antenna, for which a response can be computed,
 * but without any substructure like a beamformer
 *
 */
class Element : public Antenna {
 public:
  typedef std::shared_ptr<Element> Ptr;

  /**
   * @brief Construct a new Element object
   *
   * @param coordinate_system (antenna) CoordinateSystem
   * @param element_response ElementResponseModel
   * @param id
   */
  Element(const CoordinateSystem &coordinate_system,
          ElementResponse::Ptr element_response, int id)
      : Antenna(coordinate_system),
        id_(id),
        element_response_(element_response) {}

  Antenna::Ptr Clone() const override;

  /**
   * @brief Compute the local response of the element.
   *
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency of the plane wave (Hz).
   * @param direction Direction of arrival (ITRF, m).
   * @param id ID of element
   * @param options
   * @return matrix22c_t
   */
  virtual matrix22c_t LocalResponse(real_t time, real_t freq,
                                    const vector3r_t &direction, size_t id,
                                    const Options &options) const;

 protected:
  virtual matrix22c_t LocalResponse(real_t time, real_t freq,
                                    const vector3r_t &direction,
                                    const Options &options) const override {
    return LocalResponse(time, freq, direction, id_, options);
  };

  int id_;
  ElementResponse::Ptr element_response_;
};
}  // namespace everybeam

#endif
