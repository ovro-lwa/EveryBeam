// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_COORDS_COORDUTILS_H_
#define EVERYBEAM_COORDS_COORDUTILS_H_

#include "./../common/types.h"
#include <casacore/measures/Measures/MCDirection.h>

namespace everybeam {
namespace coords {

struct CoordinateSystem {
  std::size_t width, height;
  double ra, dec, dl, dm, phase_centre_dl, phase_centre_dm;
};

/**
 * @brief Convert Casacore itrfDir to vector3r_t
 *
 * @param itrf_dir
 * @param itrf
 */
inline void SetITRFVector(const casacore::MDirection& itrf_dir,
                          vector3r_t& itrf) {
  std::copy_n(itrf_dir.getValue().getValue().data(), itrf.size(), itrf.begin());
}
}  // namespace coords
}  // namespace everybeam
#endif  // EVERYBEAM_COORDS_COORDUTILS_H_
