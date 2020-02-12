#include "Antenna.h"

#include "MathUtil.h"

namespace LOFAR {
namespace StationResponse {

vector3r_t Antenna::transform_to_local_direction(const vector3r_t &direction) {
    vector3r_t local_direction {
        dot(m_coordinate_system.axes.p, direction),
        dot(m_coordinate_system.axes.q, direction),
        dot(m_coordinate_system.axes.r, direction),
    };

    return local_direction;
}

} // namespace StationResponse
} // namespace LOFAR