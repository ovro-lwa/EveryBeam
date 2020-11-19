// Station.cc: Representation of the station beam former.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "station.h"
#include "common/mathutils.h"
#include "beamformerlofar.h"

#include "hamaker/hamakerelementresponse.h"
#include "oskar/oskarelementresponse.h"
#include "lobes/lobeselementresponse.h"

using namespace everybeam;
using everybeam::coords::ITRFDirection;

Station::Station(const std::string &name, const vector3r_t &position,
                 const ElementResponseModel model)
    : name_(name),
      position_(position),
      phase_reference_(position),
      element_response_model_(model),
      element_response_(nullptr) {
  vector3r_t ncp = {{0.0, 0.0, 1.0}};
  ncp_.reset(new coords::ITRFDirection(ncp));
  vector3r_t ncppol0 = {{1.0, 0.0, 0.0}};
  ncp_pol0_.reset(new coords::ITRFDirection(ncppol0));

  // LOBES is currently only supported for CS302LBA. Check this.
  if (model == ElementResponseModel::kLOBES && name != "CS302LBA") {
    std::cout << "Switching to Hamaker, LOBES is not yet "
                 "supported for station "
              << name << std::endl;
    element_response_model_ = ElementResponseModel::kHamaker;
  }
  SetResponseModel(element_response_model_);
}

void Station::SetResponseModel(const ElementResponseModel model) {
  switch (model) {
    case kHamaker:
      element_response_.set(HamakerElementResponse::GetInstance(name_));
      break;
    case kOSKARDipole:
      element_response_.set(OSKARElementResponseDipole::GetInstance());
      break;
    case kOSKARSphericalWave:
      element_response_.set(OSKARElementResponseSphericalWave::GetInstance());
      break;
    case kLOBES:
      element_response_.set(LOBESElementResponse::GetInstance(name_));
      break;
    default:
      std::stringstream message;
      message << "The requested element response model '" << model
              << "' is not implemented.";
      throw std::runtime_error(message.str());
  }
}

void Station::SetResponse(std::shared_ptr<ElementResponse> element_response) {
  element_response_.set(element_response);
}

const std::string &Station::GetName() const { return name_; }

const vector3r_t &Station::GetPosition() const { return position_; }

void Station::SetPhaseReference(const vector3r_t &reference) {
  phase_reference_ = reference;
}

const vector3r_t &Station::GetPhaseReference() const {
  return phase_reference_;
}

void Station::SetAntenna(Antenna::Ptr antenna) {
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
  if (auto beamformer_lofar =
          std::dynamic_pointer_cast<BeamFormerLofar>(antenna)) {
    antenna = beamformer_lofar->GetElement();
    antenna->Transform(beamformer_lofar->coordinate_system_);
  }
  element_ = std::dynamic_pointer_cast<Element>(antenna);
}

// ========================================================
matrix22c_t Station::ComputeElementResponse(real_t time, real_t freq,
                                            const vector3r_t &direction,
                                            size_t id, bool is_local,
                                            bool rotate) const {
  Antenna::Options options;
  options.rotate = rotate;

  if (rotate) {
    vector3r_t ncp_ = NCP(time);
    vector3r_t east = normalize(cross(ncp_, direction));
    vector3r_t north = cross(direction, east);
    options.east = east;
    options.north = north;
  }

  return is_local ? element_->LocalResponse(time, freq, direction, id, options)
                  : element_->ResponseID(time, freq, direction, id, options);
}

matrix22c_t Station::ComputeElementResponse(real_t time, real_t freq,
                                            const vector3r_t &direction,
                                            bool is_local, bool rotate) const {
  Antenna::Options options;
  options.rotate = rotate;

  if (options.rotate) {
    vector3r_t ncp_ = NCP(time);
    vector3r_t east = normalize(cross(ncp_, direction));
    vector3r_t north = cross(direction, east);
    options.east = east;
    options.north = north;
  }

  return is_local ? element_->LocalResponse(time, freq, direction,
                                            element_->GetElementID(), options)
                  : element_->Response(time, freq, direction, options);
}

matrix22c_t Station::Response(real_t time, real_t freq,
                              const vector3r_t &direction, real_t freq0,
                              const vector3r_t &station0,
                              const vector3r_t &tile0,
                              const bool rotate) const {
  Antenna::Options options;
  options.freq0 = freq0;
  options.station0 = station0;
  options.tile0 = tile0;
  options.rotate = rotate;

  if (rotate) {
    vector3r_t ncp_ = NCP(time);
    vector3r_t east = normalize(cross(ncp_, direction));
    vector3r_t north = cross(direction, east);
    options.east = east;
    options.north = north;
  }

  matrix22c_t response = antenna_->Response(time, freq, direction, options);

  return response;
}

diag22c_t Station::ArrayFactor(real_t time, real_t freq,
                               const vector3r_t &direction, real_t freq0,
                               const vector3r_t &station0,
                               const vector3r_t &tile0) const {
  Antenna::Options options;
  options.freq0 = freq0;
  options.station0 = station0;
  options.tile0 = tile0;
  return antenna_->ArrayFactor(time, freq, direction, options);
}

matrix22r_t Station::Rotation(real_t time, const vector3r_t &direction) const {
  // rotation needs to be optional, normally you only want to rotate your
  // coordinatesytem for the center of your (mosaiced) image
  // Compute the cross product of the NCP and the target direction. This
  // yields a vector tangent to the celestial sphere at the target
  // direction, pointing towards the East (the direction of +Y in the IAU
  // definition, or positive right ascension).
  // Test if the direction is equal to the NCP. If it is, take a random
  // vector orthogonal to v1 (the east is not defined here).
  vector3r_t v1;
  if (std::abs(NCP(time)[0] - direction[0]) < 1e-9 &&
      std::abs(NCP(time)[1] - direction[1]) < 1e-9 &&
      std::abs(NCP(time)[2] - direction[2]) < 1e-9) {
    // Make sure v1 is orthogonal to NCP(time). In the direction of the meridian
    v1 = normalize(NCPPol0(time));
  } else {
    v1 = normalize(cross(NCP(time), direction));
  }

  // Compute the cross product of the antenna field normal (R) and the
  // target direction. This yields a vector tangent to the topocentric
  // spherical coordinate system at the target direction, pointing towards
  // the direction of positive phi (which runs East over North around the
  // pseudo zenith).
  // Test if the normal is equal to the target direction. If it is, take
  // a random vector orthogonal to the normal.
  vector3r_t v2;
  if (std::abs(antenna_->coordinate_system_.axes.r[0] - direction[0]) < 1e-9 &&
      std::abs(antenna_->coordinate_system_.axes.r[1] - direction[1]) < 1e-9 &&
      std::abs(antenna_->coordinate_system_.axes.r[2] - direction[2]) < 1e-9) {
    // Nothing to be rotated if the direction is equal to zenith
    v2 = v1;
  } else {
    v2 = normalize(cross(antenna_->coordinate_system_.axes.r, direction));
  }

  // Compute the cosine and sine of the parallactic angle, i.e. the angle
  // between v1 and v2, both tangent to a latitude circle of their
  // respective spherical coordinate systems.
  real_t coschi = dot(v1, v2);
  real_t sinchi;
  if (coschi == 1.0)
    sinchi = 0.0;
  else
    sinchi = dot(cross(v1, v2), direction);

  // The input coordinate system is a right handed system with its third
  // axis along the direction of propagation (IAU +Z). The output
  // coordinate system is right handed as well, but its third axis points
  // in the direction of arrival (i.e. exactly opposite).
  //
  // Because the electromagnetic field is always perpendicular to the
  // direction of propagation, we only need to relate the (X, Y) axes of
  // the input system to the corresponding (theta, phi) axes of the output
  // system.
  //
  // To this end, we first rotate the input system around its third axis
  // to align the Y axis with the phi axis. The X and theta axis are
  // parallel after this rotation, but point in opposite directions. To
  // align the X axis with the theta axis, we flip it.
  //
  // The Jones matrix to align the Y axis with the phi axis when these are
  // separated by an angle phi (measured counter-clockwise around the
  // direction of propagation, looking towards the origin), is given by:
  //
  // [ cos(phi)  sin(phi)]
  // [-sin(phi)  cos(phi)]
  //
  // Here, cos(phi) and sin(phi) can be computed directly, without having
  // to compute phi first (see the computation of coschi and sinchi
  // above).
  //
  // Now, sinchi as computed above is opposite to sin(phi), because the
  // direction used in the computation is the direction of arrival instead
  // of the direction of propagation. Therefore, the sign of sinchi needs
  // to be reversed. Furthermore, as explained above, the X axis has to be
  // flipped to align with the theta axis. The Jones matrix returned from
  // this function is therefore given by:
  //
  // [-coschi  sinchi]
  // [ sinchi  coschi]
  matrix22r_t rotation = {{{{-coschi, sinchi}}, {{sinchi, coschi}}}};
  return rotation;
}

vector3r_t Station::NCP(real_t time) const { return ncp_->at(time); }

vector3r_t Station::NCPPol0(real_t time) const { return ncp_pol0_->at(time); }
