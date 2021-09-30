// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "dishpoint.h"
#include "../telescope/dish.h"
#include "../circularsymmetric/voltagepattern.h"
#include "../circularsymmetric/vlabeam.h"

#include <aocommon/uvector.h>

using aocommon::UVector;

namespace everybeam {
namespace pointresponse {
void DishPoint::FullResponse(std::complex<float>* buffer, double ra, double dec,
                             double freq, [[maybe_unused]] size_t station_idx,
                             size_t field_id) {
  const telescope::Dish& dishtelescope =
      static_cast<const telescope::Dish&>(*telescope_);

  double pdir_ra = dishtelescope.ms_properties_.field_pointing[field_id].first,
         pdir_dec =
             dishtelescope.ms_properties_.field_pointing[field_id].second;
  std::array<double, 5> coefs =
      circularsymmetric::VLABeam::GetCoefficients("", freq);
  const double max_radius_arc_min = 53.0;
  circularsymmetric::VoltagePattern vp(freq, max_radius_arc_min);
  aocommon::UVector<double> coefs_vec(coefs.begin(), coefs.end());
  vp.EvaluatePolynomial(coefs_vec, false);
  vp.Render(buffer, ra, dec, pdir_ra, pdir_dec, freq);
}

void DishPoint::FullResponseAllStations(std::complex<float>* buffer, double ra,
                                        double dec, double freq,
                                        size_t field_id) {
  FullResponse(buffer, ra, dec, freq, 0., field_id);

  // Just repeat nstations times
  for (size_t i = 1; i != telescope_->GetNrStations(); ++i) {
    std::copy_n(buffer, 4, buffer + i * 4);
  }
}
}  // namespace pointresponse
}  // namespace everybeam