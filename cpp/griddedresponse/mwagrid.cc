#include "mwagrid.h"
#include "../telescope/mwa.h"

#include <aocommon/imagecoordinates.h>

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MDirection.h>

#include <memory>

using everybeam::griddedresponse::MWAGrid;
using everybeam::mwabeam::TileBeam2016;

void MWAGrid::CalculateStation(std::complex<float>* buffer, double time,
                               double frequency, size_t station_idx, size_t) {
  const telescope::MWA& mwatelescope =
      static_cast<const telescope::MWA&>(*telescope_);
  casacore::MEpoch time_epoch(casacore::Quantity(time, "s"));
  casacore::MeasFrame frame(mwatelescope.ms_properties_.array_position,
                            time_epoch);

  const casacore::MDirection::Ref hadec_ref(casacore::MDirection::HADEC, frame);
  const casacore::MDirection::Ref azelgeo_ref(casacore::MDirection::AZELGEO,
                                              frame);
  const casacore::MDirection::Ref j2000_ref(casacore::MDirection::J2000, frame);
  casacore::MDirection::Convert j2000_to_hadecref(j2000_ref, hadec_ref),
      j2000_to_azelgeoref(j2000_ref, azelgeo_ref);
  casacore::MPosition wgs = casacore::MPosition::Convert(
      mwatelescope.ms_properties_.array_position, casacore::MPosition::WGS84)();
  double arrLatitude = wgs.getValue().getLat();

  if (!tile_beam_) {
    tile_beam_.reset(
        new TileBeam2016(mwatelescope.ms_properties_.delays,
                         mwatelescope.GetOptions().frequency_interpolation,
                         mwatelescope.GetOptions().coeff_path));
  }
  std::complex<float>* buffer_ptr = buffer;
  for (size_t y = 0; y != height_; ++y) {
    for (size_t x = 0; x != width_; ++x) {
      double l, m, ra, dec;
      aocommon::ImageCoordinates::XYToLM(x, y, dl_, dm_, width_, height_, l, m);
      l += phase_centre_dl_;
      m += phase_centre_dm_;
      aocommon::ImageCoordinates::LMToRaDec(l, m, ra_, dec_, ra, dec);

      std::complex<double> gain[4];
      tile_beam_->ArrayResponse(ra, dec, j2000_ref, j2000_to_hadecref,
                                j2000_to_azelgeoref, arrLatitude, frequency,
                                gain);

      for (size_t i = 0; i != 4; ++i) {
        *buffer_ptr = gain[i];
        ++buffer_ptr;
      }
    }
  }
};

void MWAGrid::CalculateAllStations(std::complex<float>* buffer, double time,
                                   double frequency, size_t) {
  CalculateStation(buffer, time, frequency, 0, 0);
  // Repeated copy for nstations
  for (size_t i = 1; i != telescope_->GetNrStations(); ++i) {
    std::copy_n(buffer, width_ * height_ * 4,
                buffer + i * width_ * height_ * 4);
  }
};