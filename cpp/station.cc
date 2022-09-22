// Station.cc: Representation of the station beam former.
//
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "station.h"
#include "common/mathutils.h"
#include "beamformerlofar.h"

using namespace everybeam;
using everybeam::coords::ITRFDirection;

Station::Station(const std::string& name, const vector3r_t& position,
                 const Options& options)
    : name_(name),
      position_(position),
      options_(options),
      phase_reference_(position),
      element_response_(ElementResponse::GetInstance(
          options.element_response_model, name_, options_)) {
  const vector3r_t ncp = {0.0, 0.0, 1.0};
  ncp_.reset(new coords::ITRFDirection(ncp));
  const vector3r_t ncppol0 = {1.0, 0.0, 0.0};
  ncp_pol0_.reset(new coords::ITRFDirection(ncppol0));
}

const std::string& Station::GetName() const { return name_; }

const vector3r_t& Station::GetPosition() const { return position_; }

void Station::SetPhaseReference(const vector3r_t& reference) {
  phase_reference_ = reference;
}

const vector3r_t& Station::GetPhaseReference() const {
  return phase_reference_;
}

void Station::SetAntenna(std::shared_ptr<Antenna> antenna) {
  antenna_ = antenna;

  // The antenna can be either an Element or a BeamFormer
  // If it is a BeamFormer we recursively extract the first antenna
  // until we have a BeamFormerLofar or an Element.
  //
  // The extraction returns copies so antenna_ remains unchanged.
  // The element that is found is used in ComputeElementResponse to
  // compute the element response.

  while (auto beamformer = std::dynamic_pointer_cast<BeamFormer>(antenna)) {
    antenna = beamformer->ExtractAntenna(0);
  }

  // If we can cast to BeamFormerLofar, then extract the Element - please
  // note that the Element was upcasted from an ElementHamaker into an Element
  // in BeamFormerLofarHBA/LBA::Clone()!- and Transform the Element with the
  // coordinate system of the HBA/LBA beam former.
  if (auto beamformer_lofar = dynamic_cast<BeamFormerLofar*>(antenna.get())) {
    element_ = beamformer_lofar->GetElement();
    element_->Transform(beamformer_lofar->coordinate_system_);
  } else {
    element_ = std::dynamic_pointer_cast<Element>(antenna);
  }
}

// ========================================================
aocommon::MC2x2 Station::ComputeElementResponse(real_t time, real_t freq,
                                                const vector3r_t& direction,
                                                size_t id, bool is_local,
                                                bool rotate) const {
  Antenna::Options options;
  options.rotate = rotate;

  if (rotate) {
    const vector3r_t ncp_t = NCP(time);
    const vector3r_t east = normalize(cross(ncp_t, direction));
    const vector3r_t north = cross(direction, east);
    options.east = east;
    options.north = north;
  }

  return is_local ? element_->LocalResponse(*element_response_, time, freq,
                                            direction, id, options)
                  : element_->ResponseID(*element_response_, time, freq,
                                         direction, id, options);
}

aocommon::MC2x2 Station::ComputeElementResponse(real_t time, real_t freq,
                                                const vector3r_t& direction,
                                                bool is_local,
                                                bool rotate) const {
  return ComputeElementResponse(time, freq, direction, element_->GetElementID(),
                                is_local, rotate);
}

aocommon::MC2x2 Station::Response(real_t time, real_t freq,
                                  const vector3r_t& direction, real_t freq0,
                                  const vector3r_t& station0,
                                  const vector3r_t& tile0,
                                  const bool rotate) const {
  Antenna::Options options;
  options.freq0 = freq0;
  options.station0 = station0;
  options.tile0 = tile0;
  options.rotate = rotate;

  if (rotate) {
    const vector3r_t ncp_t = NCP(time);
    const vector3r_t east = normalize(cross(ncp_t, direction));
    const vector3r_t north = cross(direction, east);
    options.east = east;
    options.north = north;
  }

  return antenna_->Response(*element_response_, time, freq, direction, options);
}

aocommon::MC2x2Diag Station::ArrayFactor(real_t time, real_t freq,
                                         const vector3r_t& direction,
                                         real_t freq0,
                                         const vector3r_t& station0,
                                         const vector3r_t& tile0) const {
  Antenna::Options options;
  options.freq0 = freq0;
  options.station0 = station0;
  options.tile0 = tile0;
  return antenna_->ArrayFactor(time, freq, direction, options);
}

vector3r_t Station::NCP(real_t time) const { return ncp_->at(time); }

vector3r_t Station::NCPPol0(real_t time) const { return ncp_pol0_->at(time); }
