#ifndef ANTENNA_H
#define ANTENNA_H

#include <complex>
#include <memory>

#include "Types.h"

namespace LOFAR {
namespace StationResponse {

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
    };



    typedef std::shared_ptr<Antenna> Ptr;

    struct Options
    {
        real_t freq0;
        const vector3r_t *station0;
        const vector3r_t *tile0;
    };

    Antenna() {}

    Antenna(CoordinateSystem &coordinate_system) :
        m_coordinate_system(coordinate_system)
    {
    }

    matrix22c_t response(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options = {})
    {
        // TODO
        // Transform direction and directions in options to local coordinatesystem
        return local_response(time, freq, direction, options);
    }

    diag22c_t arrayFactor(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options = {})
    {
        // TODO
        // Transform direction and directions in options to local coordinatesystem
        return local_arrayFactor(time, freq, direction, options);
    }

    CoordinateSystem m_coordinate_system;
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

};

} // namespace StationResponse
} // namespace LOFAR
#endif
