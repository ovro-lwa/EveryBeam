#ifndef EVERYBEAM_COORDS_COORDUTILS_H_
#define EVERYBEAM_COORDS_COORDUTILS_H_

#include "./../common/Types.h"
#include <casacore/measures/Measures/MCDirection.h>

namespace everybeam {
namespace coords {

/**
 * @brief Convert Casacore itrfDir to vector3r_t
 *
 * @param itrfDir
 * @param itrf
 */
void setITRFVector(const casacore::MDirection& itrfDir, vector3r_t& itrf) {
  const casacore::Vector<double>& itrfVal = itrfDir.getValue().getValue();
  itrf[0] = itrfVal[0];
  itrf[1] = itrfVal[1];
  itrf[2] = itrfVal[2];
}
}  // namespace coords
}  // namespace everybeam
#endif  // EVERYBEAM_COORDS_COORDUTILS_H_