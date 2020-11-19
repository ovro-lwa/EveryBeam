// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mwagrid.h"
#include "../telescope/mwa.h"

#include <aocommon/imagecoordinates.h>

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MDirection.h>

#include <memory>

using aocommon::HMC4x4;
using aocommon::MC2x2;
using aocommon::UVector;
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
}

void MWAGrid::CalculateAllStations(std::complex<float>* buffer, double time,
                                   double frequency, size_t) {
  CalculateStation(buffer, time, frequency, 0, 0);
  // Repeated copy for nstations
  for (size_t i = 1; i != telescope_->GetNrStations(); ++i) {
    std::copy_n(buffer, width_ * height_ * 4,
                buffer + i * width_ * height_ * 4);
  }
}

void MWAGrid::MakeIntegratedSnapshot(std::vector<aocommon::HMC4x4>& matrices,
                                     double time, double frequency,
                                     size_t field_id,
                                     const double* baseline_weights_interval) {
  size_t nstations = telescope_->GetNrStations();
  UVector<std::complex<float>> buffer_undersampled(
      GetStationBufferSize(nstations));
  CalculateAllStations(buffer_undersampled.data(), time, frequency, field_id);

  // For MWA, we can simply weight a (time) snapshot with the accumulated
  // baseline weights
  size_t nbaselines = nstations * (nstations + 1) / 2;
  double snapshot_weight = 0.;
  for (size_t index = 0; index != nbaselines; ++index) {
    snapshot_weight += baseline_weights_interval[index];
  }

  for (size_t y = 0; y != height_; ++y) {
    for (size_t x = 0; x != width_; ++x) {
      size_t offset = (y * width_ + x) * 4;
      MC2x2 A(buffer_undersampled[offset], buffer_undersampled[offset + 1],
              buffer_undersampled[offset + 2], buffer_undersampled[offset + 3]);

      // Mueller matrix constant for all baselines, so just compute once for
      // each individual pixel
      matrices[y * width_ + x] =
          HMC4x4::KroneckerProduct(A.HermTranspose().Transpose(), A) *
          snapshot_weight;
    }
  }
}