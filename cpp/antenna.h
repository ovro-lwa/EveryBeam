// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ANTENNA_H
#define EVERYBEAM_ANTENNA_H

#include <complex>
#include <memory>
#include <iostream>

#include <aocommon/matrix2x2.h>
#include <aocommon/matrix2x2diag.h>

#include "common/types.h"

namespace everybeam {

/**
 * @brief (Virtual) class describing an antenna, and computing the corresponding
 * Response() and ArrayFactor().
 *
 */
class Antenna {
 public:
  /**
   *  \brief %Station coordinate system.
   *
   *  A right handed, cartesian, local coordinate system with coordinate axes
   *  \p p, \p q, and \p r is associated with each antenna field.
   *
   *  The r-axis is orthogonal to the antenna field, and points towards the
   *  local pseudo zenith.
   *
   *  The q-axis is the northern bisector of the \p X and \p Y dipoles, i.e.
   *  it is the reference direction from which the orientation of the dual
   *  dipole antennae is determined. The q-axis points towards the North at
   *  the core. At remote sites it is defined as the intersection of the
   *  antenna field plane and a plane parallel to the meridian plane at the
   *  core. This ensures the reference directions at all sites are similar.
   *
   *  The p-axis is orthogonal to both other axes, and points towards the East
   *  at the core.
   *
   *  The axes and origin of the anntena field coordinate system are expressed
   *  as vectors in the geocentric, cartesian, ITRF coordinate system, in
   *  meters.
   *
   *  \sa "LOFAR Reference Plane and Reference Direction", M.A. Brentjens,
   *  LOFAR-ASTRON-MEM-248.
   */
  struct CoordinateSystem {
    struct Axes {
      vector3r_t p;
      vector3r_t q;
      vector3r_t r;
    };
    vector3r_t origin;
    Axes axes;

    constexpr static Axes identity_axes =
        Axes{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    constexpr static vector3r_t zero_origin = vector3r_t{0.0, 0.0, 0.0};
  };

  constexpr static CoordinateSystem IdentityCoordinateSystem{
      CoordinateSystem::zero_origin, CoordinateSystem::identity_axes};

  typedef std::shared_ptr<Antenna> Ptr;

  /**
   * @brief Struct containing antenna options
   *
   */
  struct Options {
    real_t freq0;         //!< %Antenna reference frequency (Hz).
    vector3r_t station0;  //!< Reference direction (ITRF, m)
    vector3r_t tile0;     //!< Tile beam former reference direction (ITRF, m).
    bool
        rotate;  //!< Boolean deciding if paralactic rotation should be applied.
    vector3r_t east;   //!< Eastward pointing unit vector
    vector3r_t north;  //!< Northward pointing unit vector
  };

  /**
   * @brief Construct a new %Antenna object
   *
   */
  Antenna()
      :  // default coordinate system
         // no shift of origin, no rotation
        Antenna({CoordinateSystem::zero_origin,  // origin
                 CoordinateSystem::identity_axes}) {}

  /**
   * @brief Construct a new %Antenna object, given a coordinate system
   *
   * @param coordinate_system
   */
  Antenna(const CoordinateSystem &coordinate_system)
      :  // default phase reference system is the origin of the coordinate
         // system
        Antenna(coordinate_system, coordinate_system.origin) {}

  virtual ~Antenna(){};

  /**
   * @brief Construct a new %Antenna object, given a coordinate system and a
   * phase reference position.
   *
   * @param coordinate_system Coordinate system
   * @param phase_reference_position Phase reference position
   */
  Antenna(const CoordinateSystem &coordinate_system,
          const vector3r_t &phase_reference_position)
      : coordinate_system_(coordinate_system),
        phase_reference_position_(phase_reference_position),
        enabled_{true, true} {}

  Antenna(const vector3r_t &phase_reference_position)
      : coordinate_system_({phase_reference_position,  // origin
                            CoordinateSystem::identity_axes}),
        phase_reference_position_(phase_reference_position),
        enabled_{true, true} {}

  /**
   * @brief Makes a copy of this Antenna object
   *
   * The method is virtual, so that copies can be created from a pointer
   * to the base (Antenna) class.
   * The original remains unchanged, therefore the method is const.
   * The method has no implementation in the Antenna class, because
   * Antenna is abstract, so no copy can be instantiated.
   *
   * This method is used by the ExtractAntenna method of the BeamFormer
   * class to create a copy of one of the Antennas it contains.
   */
  virtual Ptr Clone() const = 0;

  /**
   * @brief Transform internal coordinate systems and positions
   *
   * @param coordinate_system to apply in the transformation
   *
   * This method is used by BeamFormer::ExtractAntenna to lift
   * an antenna out of the beamformer.
   *
   * The transformation is needed because the coordinate system of
   * an antenna in a beamformer is expressed in terms of
   * the coordinate system of the beamformer.
   * To turn an embedded antenna into a stand-alone antenna,
   * the coordinate system of the beamformer needs to be
   * applied to the coordinate system of the antenna
   */
  void Transform(const CoordinateSystem &coordinate_system);

  /**
   * @brief Compute the %Antenna Response
   *
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency of the plane wave (Hz).
   * @param direction Direction of arrival (ITRF, m).
   * @param options
   */
  virtual aocommon::MC2x2 Response(real_t time, real_t freq,
                                   const vector3r_t &direction,
                                   const Options &options = {}) const {
    // Transform direction and directions in options to local coordinatesystem
    vector3r_t local_direction = TransformToLocalDirection(direction);
    Options local_options;
    local_options.freq0 = options.freq0;
    local_options.station0 = TransformToLocalDirection(options.station0);
    local_options.tile0 = TransformToLocalDirection(options.tile0);
    local_options.rotate = options.rotate;
    local_options.east = TransformToLocalDirection(options.east);
    local_options.north = TransformToLocalDirection(options.north);
    return LocalResponse(time, freq, local_direction, local_options);
  }

  /**
   * @brief Compute the array factor of the antenna
   *
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency of the plane wave (Hz).
   * @param direction Direction of arrival (ITRF, m).
   * @param options
   */
  virtual aocommon::MC2x2Diag ArrayFactor(real_t time, real_t freq,
                                          const vector3r_t &direction,
                                          const Options &options) const {
    // Transform direction and directions in options to local coordinatesystem
    const vector3r_t local_direction = TransformToLocalDirection(direction);
    Options local_options;
    local_options.freq0 = options.freq0;
    local_options.station0 = TransformToLocalDirection(options.station0);
    local_options.tile0 = TransformToLocalDirection(options.tile0);
    return LocalArrayFactor(time, freq, local_direction, local_options);
  }

  CoordinateSystem coordinate_system_;
  vector3r_t phase_reference_position_;
  bool enabled_[2];

 protected:
  vector3r_t TransformToLocalDirection(const vector3r_t &direction) const;

 private:
  virtual aocommon::MC2x2 LocalResponse(real_t time, real_t freq,
                                        const vector3r_t &direction,
                                        const Options &options) const = 0;

  virtual aocommon::MC2x2Diag LocalArrayFactor(
      [[maybe_unused]] real_t time, [[maybe_unused]] real_t freq,
      [[maybe_unused]] const vector3r_t &direction,
      [[maybe_unused]] const Options &options) const {
    return aocommon::MC2x2Diag::Unity();
  }
};

}  // namespace everybeam
#endif
