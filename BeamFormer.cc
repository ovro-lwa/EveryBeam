#include "BeamFormer.h"

#include <cmath>

constexpr double speed_of_light = 299792458.0;

namespace LOFAR {
namespace StationResponse {

vector3r_t BeamFormer::transform_to_local_position(const vector3r_t &position) {
    vector3r_t dposition {
        position[0] - m_coordinate_system.origin[0],
        position[1] - m_coordinate_system.origin[1],
        position[2] - m_coordinate_system.origin[2]
    };
    vector3r_t local_position {
        dot(m_coordinate_system.axes.p, dposition),
        dot(m_coordinate_system.axes.q, dposition),
        dot(m_coordinate_system.axes.r, dposition),
    };

    return local_position;
}


std::vector<std::complex<double>> BeamFormer::compute_geometric_response(double freq, const vector3r_t &direction) const
{
    std::vector<std::complex<double>> result;
    result.reserve(m_antennas.size());
    for (auto &antenna : m_antennas)
    {
        double dl = direction[0] * antenna->m_coordinate_system.origin[0] +
                    direction[1] * antenna->m_coordinate_system.origin[1] +
                    direction[2] * antenna->m_coordinate_system.origin[2];

        double phase = 2 * M_PI * dl * speed_of_light * freq;
        result.push_back({std::sin(phase), std::cos(phase)});
    }
    return result;
}

vector3r_t BeamFormer::compute_local_pointing(double time) const
{
    return {m_pointing[0], m_pointing[1], m_pointing[1]};
}


std::vector<std::pair<std::complex<double>,std::complex<double>>> BeamFormer::compute_weights(double time, double freq) const
{
    std::vector<std::pair<std::complex<double>,std::complex<double>>> result;

    vector3r_t pointing = compute_local_pointing(time);

    auto geometric_response = compute_geometric_response(freq, pointing);
    result.reserve(geometric_response.size());
    for (auto phasor : geometric_response)
    {
        result.push_back({std::conj(phasor), std::conj(phasor)});
    }
    return result;
}


matrix22c_t BeamFormer::local_response(
    real_t time,
    real_t freq,
    const vector3r_t &direction,
    const Options &options) const
{
    auto weights = compute_weights(time, freq);
    auto geometric_response = compute_geometric_response(freq, direction);

    matrix22c_t result = {0};
    for (unsigned int antenna_idx = 0; antenna_idx < m_antennas.size(); ++antenna_idx) {
        auto antenna = m_antennas[antenna_idx];
        auto antenna_weight = weights[antenna_idx];
        auto antenna_geometric_reponse = geometric_response[antenna_idx];

        matrix22c_t antenna_response = antenna->response(time, freq, direction, options);

        result[0][0] += antenna_weight.first * antenna_geometric_reponse * antenna_response[0][0];
        result[0][1] += antenna_weight.first * antenna_geometric_reponse * antenna_response[0][1];
        result[1][0] += antenna_weight.second * antenna_geometric_reponse * antenna_response[1][0];
        result[1][1] += antenna_weight.second * antenna_geometric_reponse * antenna_response[1][1];

    }
    return result;
}


} // namespace StationResponse
} // namespace LOFAR

