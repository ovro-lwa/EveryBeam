#ifndef EVERYBEAM_ELEMENT_H
#define EVERYBEAM_ELEMENT_H

#include <complex>
#include <memory>

#include "Antenna.h"
#include "ElementResponse.h"
#include "common/Types.h"

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
        m_id(id),
        m_element_response(element_response) {}

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
  matrix22c_t local_response(real_t time, real_t freq,
                             const vector3r_t &direction, size_t id,
                             const Options &options) const;

 private:
  virtual matrix22c_t local_response(
      real_t time, real_t freq, const vector3r_t &direction,
      const Options &options) const final override;

  int m_id;
  ElementResponse::Ptr m_element_response;
};

}  // namespace everybeam

#endif
