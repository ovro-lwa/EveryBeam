// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "dishpoint.h"
#include "../telescope/dish.h"
#include "../circularsymmetric/voltagepattern.h"
#include "../circularsymmetric/vlacoefficients.h"

#include <aocommon/uvector.h>

namespace everybeam {
namespace pointresponse {
void DishPoint::Response(BeamMode /* beam_mode */, std::complex<float>* buffer,
                         double ra, double dec, double freq,
                         size_t /* station_idx */, size_t field_id) {
  const telescope::Dish& dish_telescope =
      static_cast<const telescope::Dish&>(*telescope_);

  double pdir_ra;
  double pdir_dec;
  std::tie(pdir_ra, pdir_dec) = dish_telescope.GetFieldPointing()[field_id];
  const double max_radius_arc_min =
      dish_telescope.GetDishCoefficients()->MaxRadiusInArcMin();
  const double reference_frequency =
      dish_telescope.GetDishCoefficients()->ReferenceFrequency();
  circularsymmetric::VoltagePattern vp(
      dish_telescope.GetDishCoefficients()->GetFrequencies(freq),
      max_radius_arc_min, reference_frequency);
  const aocommon::UVector<double> coefs_vec =
      dish_telescope.GetDishCoefficients()->GetCoefficients(freq);
  vp.EvaluatePolynomial(coefs_vec,
                        dish_telescope.GetDishCoefficients()->AreInverted());
  vp.Render(buffer, ra, dec, pdir_ra, pdir_dec, freq);
}

void DishPoint::ResponseAllStations(BeamMode beam_mode,
                                    std::complex<float>* buffer, double ra,
                                    double dec, double freq, size_t field_id) {
  Response(beam_mode, buffer, ra, dec, freq, 0u, field_id);

  // Just repeat nstations times
  for (size_t i = 1; i != telescope_->GetNrStations(); ++i) {
    std::copy_n(buffer, 4, buffer + i * 4);
  }
}
}  // namespace pointresponse
}  // namespace everybeam
