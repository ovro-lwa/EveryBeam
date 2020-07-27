#include "dishgrid.h"
#include "../telescope/dish.h"
#include "../circularsymmetric/voltagepattern.h"
#include "../circularsymmetric/vlabeam.h"

#include <aocommon/uvector.h>

#include <algorithm>

using namespace everybeam;
using namespace everybeam::griddedresponse;

void DishGrid::CalculateStation(std::complex<float>* buffer, double,
                                double frequency, const size_t,
                                const size_t field_id) {
  const telescope::Dish& dishtelescope =
      static_cast<const telescope::Dish&>(*telescope_);

  double pdir_ra = dishtelescope.ms_properties_.field_pointing[field_id].first,
         pdir_dec =
             dishtelescope.ms_properties_.field_pointing[field_id].second;
  std::array<double, 5> coefs =
      circularsymmetric::VLABeam::GetCoefficients("", frequency);
  circularsymmetric::VoltagePattern vp(frequency, 53.0);
  aocommon::UVector<double> coefs_vec(coefs.begin(), coefs.end());
  vp.EvaluatePolynomial(coefs_vec, false);
  vp.Render(buffer, width_, height_, dl_, dm_, ra_, dec_, pdir_ra, pdir_dec,
            phase_centre_dl_, phase_centre_dm_, frequency);
};

void DishGrid::CalculateAllStations(std::complex<float>* buffer, double,
                                    double frequency, const size_t field_id) {
  CalculateStation(buffer, 0., frequency, 0, field_id);

  // Just repeat nstations times
  for (size_t i = 1; i != telescope_->GetNrStations(); ++i) {
    std::copy_n(buffer, width_ * height_ * 4,
                buffer + i * width_ * height_ * 4);
  }
}
