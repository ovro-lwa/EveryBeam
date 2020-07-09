// Dir2ITRF.h: Convertor that maps time to an ITRF direction.

#ifndef EVERYBEAM_DIR2ITRF_H
#define EVERYBEAM_DIR2ITRF_H

// \file
// Functor that maps J2000 to an ITRF direction.

#include "./../common/types.h"

#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCDirection.h>

#include <memory>

namespace everybeam {
namespace coords {
/**
 * @brief Class providing utilities for coordinate transformations
 * to and from ITRF (International Terrestrial Reference Frame)
 *
 */
class ITRFConverter {
 public:
  typedef std::unique_ptr<ITRFDirection> Ptr;
  typedef std::unique_ptr<const ITRFDirection> ConstPtr;

  ITRFConverter(real_t time);

  void setTime(real_t time);
  vector3r_t j2000ToITRF(const vector2r_t &j2000Direction) const;
  vector3r_t j2000ToITRF(const vector3r_t &j2000Direction) const;
  vector3r_t toITRF(const casacore::MDirection &direction) const;
  casacore::MDirection toDirection(const vector2r_t &j2000Direction) const;
  casacore::MDirection toDirection(const vector3r_t &j2000Direction) const;
  casacore::MDirection toDirection(const casacore::MDirection &direction) const;

 private:
  casacore::MeasFrame itsFrame;
  mutable casacore::MDirection::Convert itsConverter;
};
}  // namespace coords
}  // namespace everybeam
#endif
