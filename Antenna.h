#ifndef EVERYBEAM_ANTENNA_H
#define EVERYBEAM_ANTENNA_H

#include <complex>
#include <memory>
#include <iostream>

#include "Types.h"

namespace everybeam {

class Antenna
{
public:

    /**
     *  \brief Station coordinate system.
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
    struct CoordinateSystem
    {
        struct Axes
        {
            vector3r_t  p;
            vector3r_t  q;
            vector3r_t  r;
        };
        vector3r_t  origin;
        Axes        axes;

        constexpr static Axes identity_axes = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };

        constexpr static vector3r_t zero_origin = {0.0, 0.0, 0.0};
    };

    constexpr static CoordinateSystem identity_coordinate_system{
        CoordinateSystem::zero_origin,
        CoordinateSystem::identity_axes
    };

    typedef std::shared_ptr<Antenna> Ptr;

    struct Options
    {
        real_t freq0;
        vector3r_t station0;
        vector3r_t tile0;
        bool rotate;
        vector3r_t east;
        vector3r_t north;
    };

    Antenna() :
        // default coordinate system
        // no shift of origin, no rotation
        Antenna({
            CoordinateSystem::zero_origin, // origin
            CoordinateSystem::identity_axes
        })
    {}

    Antenna(const CoordinateSystem &coordinate_system) :
        // default phase reference system is the origin of the coordinate system
        Antenna(coordinate_system, coordinate_system.origin)
    {}

    Antenna(const CoordinateSystem &coordinate_system, const vector3r_t &phase_reference_position) :
        m_coordinate_system(coordinate_system),
        m_phase_reference_position(phase_reference_position),
        m_enabled{true, true}
    {
    }

    matrix22c_t response(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options = {})
    {
        // Transform direction and directions in options to local coordinatesystem
        vector3r_t local_direction = transform_to_local_direction(direction);
        Options local_options = {
            .freq0 = options.freq0,
            .station0 = transform_to_local_direction(options.station0),
            .tile0 = transform_to_local_direction(options.tile0),
            .rotate = options.rotate,
            .east = transform_to_local_direction(options.east),
            .north = transform_to_local_direction(options.north)
        };
        matrix22c_t response = local_response(time, freq, local_direction, local_options);
        return response;
    }

    diag22c_t arrayFactor(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options = {})
    {
        // Transform direction and directions in options to local coordinatesystem
        vector3r_t local_direction = transform_to_local_direction(direction);
        Options local_options = {
            .freq0 = options.freq0,
            .station0 = transform_to_local_direction(options.station0),
            .tile0 = transform_to_local_direction(options.tile0)
        };
        return local_arrayFactor(time, freq, local_direction, local_options);
    }

    CoordinateSystem m_coordinate_system;
    vector3r_t  m_phase_reference_position;
    bool m_enabled[2];

private:

    virtual matrix22c_t local_response(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options) const = 0;

    virtual diag22c_t local_arrayFactor(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options) const
    { return {1.0, 1.0}; }

    vector3r_t transform_to_local_direction(const vector3r_t &direction);

};

} // namespace everybeam
#endif
