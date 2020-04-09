#include "Element.h"
#include "MathUtil.h"

namespace LOFAR {
namespace StationResponse {

matrix22c_t Element::local_response(
    real_t time,
    real_t freq,
    const vector3r_t &direction,
    size_t id,
    const Options &options) const
{
    vector2r_t thetaphi = cart2thetaphi(direction);

    matrix22c_t result;
    static_assert(sizeof(std::complex<double>[2][2]) == sizeof(matrix22c_t));
    m_element_response->response(id, freq, thetaphi[0], thetaphi[1], reinterpret_cast<std::complex<double> (&)[2][2]>(result));

    if (options.rotate) {
        vector3r_t up = {0.0, 0.0, 1.0};
        vector3r_t e_phi = normalize(cross(up, direction));
        vector3r_t e_theta = cross(e_phi, direction);
        matrix22r_t rotation;
        rotation[0] = {dot(e_theta, options.north), dot(e_theta, options.east)};
        rotation[1] = {dot(e_phi, options.north), dot(e_phi, options.east)};
        result = result * rotation;
    }
    return result;
}

matrix22c_t Element::local_response(
    real_t time,
    real_t freq,
    const vector3r_t &direction,
    const Options &options) const
{
    return local_response(time, freq, direction, m_id, options);
}

} // namespace StationResponse
} // namespace LOFAR