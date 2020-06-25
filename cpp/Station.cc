// Station.cc: Representation of the station beam former.
//
// Copyright (C) 2013
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//
// $Id$

#include "Station.h"
#include "common/MathUtil.h"

#include "hamaker/HamakerElementResponse.h"
#include "oskar/OSKARElementResponse.h"
#include "lobes/LOBESElementResponse.h"

using namespace everybeam;

Station::Station(const std::string &name, const vector3r_t &position,
                 const ElementResponseModel model)
    : itsName(name),
      itsPosition(position),
      itsPhaseReference(position),
      itsElementResponse(nullptr) {
  setModel(model);
  vector3r_t ncp = {{0.0, 0.0, 1.0}};
  itsNCP.reset(new coords::ITRFDirection(ncp));
  vector3r_t ncppol0 = {{1.0, 0.0, 0.0}};
  itsNCPPol0.reset(new coords::ITRFDirection(ncppol0));
}

void Station::setModel(const ElementResponseModel model) {
  switch (model) {
    case Hamaker:
      itsElementResponse.set(HamakerElementResponse::getInstance(itsName));
      break;
    case OSKARDipole:
      itsElementResponse.set(OSKARElementResponseDipole::getInstance());
      break;
    case OSKARSphericalWave:
      itsElementResponse.set(OSKARElementResponseSphericalWave::getInstance());
      break;
    case LOBES:
      itsElementResponse.set(LOBESElementResponse::getInstance(itsName));
      break;
    default:
      std::stringstream message;
      message << "The requested element response model '" << model
              << "' is not implemented.";
      throw std::runtime_error(message.str());
  }
}

const std::string &Station::name() const { return itsName; }

const vector3r_t &Station::position() const { return itsPosition; }

void Station::setPhaseReference(const vector3r_t &reference) {
  itsPhaseReference = reference;
}

const vector3r_t &Station::phaseReference() const { return itsPhaseReference; }

// ========================================================
matrix22c_t Station::elementResponse(real_t time, real_t freq,
                                     const vector3r_t &direction, size_t id,
                                     const bool rotate) const {
  Antenna::Options options;
  options.rotate = rotate;

  if (rotate) {
    vector3r_t ncp_ = ncp(time);
    vector3r_t east = normalize(cross(ncp_, direction));
    vector3r_t north = cross(direction, east);
    options.east = east;
    options.north = north;
  }

  return itsElement->local_response(time, freq, direction, id, options);
}

matrix22c_t Station::elementResponse(real_t time, real_t freq,
                                     const vector3r_t &direction,
                                     const bool rotate) const {
  //     if (rotate)
  //       return itsElement->response(time, freq, direction)
  //           * rotation(time, direction);
  //     else
  //       return itsElement->response(time, freq, direction);

  return itsElement->response(time, freq, direction);
}

matrix22c_t Station::response(real_t time, real_t freq,
                              const vector3r_t &direction, real_t freq0,
                              const vector3r_t &station0,
                              const vector3r_t &tile0,
                              const bool rotate) const {
  Antenna::Options options = {
      .freq0 = freq0, .station0 = station0, .tile0 = tile0, .rotate = rotate};

  if (rotate) {
    vector3r_t ncp_ = ncp(time);
    vector3r_t east = normalize(cross(ncp_, direction));
    vector3r_t north = cross(direction, east);
    options.east = east;
    options.north = north;
  }

  matrix22c_t response = itsAntenna->response(time, freq, direction, options);

  //     if (rotate) {
  //         std::cout << "rotate" << std::endl;
  //         auto r = rotation(time, direction);
  //         std::cout << r[0][0] << ", " << r[0][1] << std::endl;
  //         std::cout << r[1][0] << ", " << r[1][1] << std::endl;
  //         response = response * r;
  //         std::cout << response[0][0] << std::endl;
  //     }
  return response;
}

diag22c_t Station::arrayFactor(real_t time, real_t freq,
                               const vector3r_t &direction, real_t freq0,
                               const vector3r_t &station0,
                               const vector3r_t &tile0) const {
  Antenna::Options options = {
      .freq0 = freq0, .station0 = station0, .tile0 = tile0};
  return itsAntenna->arrayFactor(time, freq, direction, options);
}

matrix22r_t Station::rotation(real_t time, const vector3r_t &direction) const {
  // rotation needs to be optional, normally you only want to rotate your
  // coordinatesytem for the center of your (mosaiced) image
  // Compute the cross product of the NCP and the target direction. This
  // yields a vector tangent to the celestial sphere at the target
  // direction, pointing towards the East (the direction of +Y in the IAU
  // definition, or positive right ascension).
  // Test if the direction is equal to the NCP. If it is, take a random
  // vector orthogonal to v1 (the east is not defined here).
  vector3r_t v1;
  if (std::abs(ncp(time)[0] - direction[0]) < 1e-9 &&
      std::abs(ncp(time)[1] - direction[1]) < 1e-9 &&
      std::abs(ncp(time)[2] - direction[2]) < 1e-9) {
    // Make sure v1 is orthogonal to ncp(time). In the direction of the meridian
    v1 = normalize(ncppol0(time));
  } else {
    v1 = normalize(cross(ncp(time), direction));
  }

  // Compute the cross product of the antenna field normal (R) and the
  // target direction. This yields a vector tangent to the topocentric
  // spherical coordinate system at the target direction, pointing towards
  // the direction of positive phi (which runs East over North around the
  // pseudo zenith).
  // Test if the normal is equal to the target direction. If it is, take
  // a random vector orthogonal to the normal.
  vector3r_t v2;
  if (std::abs(itsAntenna->m_coordinate_system.axes.r[0] - direction[0]) <
          1e-9 &&
      std::abs(itsAntenna->m_coordinate_system.axes.r[1] - direction[1]) <
          1e-9 &&
      std::abs(itsAntenna->m_coordinate_system.axes.r[2] - direction[2]) <
          1e-9) {
    // Nothing to be rotated if the direction is equal to zenith
    v2 = v1;
  } else {
    v2 = normalize(cross(itsAntenna->m_coordinate_system.axes.r, direction));
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

vector3r_t Station::ncp(real_t time) const { return itsNCP->at(time); }

vector3r_t Station::ncppol0(real_t time) const { return itsNCPPol0->at(time); }
