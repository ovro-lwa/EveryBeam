#include "BeamFormer.h"

#include "common/Constants.h"
#include "common/MathUtil.h"

#include <cmath>

namespace everybeam {
vector3r_t BeamFormer::transform_to_local_position(const vector3r_t &position) {
  // Get antenna position relative to coordinate system origin
  vector3r_t dposition{position[0] - m_coordinate_system.origin[0],
                       position[1] - m_coordinate_system.origin[1],
                       position[2] - m_coordinate_system.origin[2]};
  // Inner product on orthogonal unit vectors of coordinate system
  vector3r_t local_position{
      dot(m_coordinate_system.axes.p, dposition),
      dot(m_coordinate_system.axes.q, dposition),
      dot(m_coordinate_system.axes.r, dposition),
  };

  return local_position;
}

std::vector<std::complex<double>> BeamFormer::compute_geometric_response(
    const double freq, const vector3r_t &direction) const {
  // Initialize and fill result vector by looping over antennas
  std::vector<std::complex<double>> result(m_antennas.size());
  for (auto &antenna : m_antennas) {
    const double dl = direction[0] * (antenna->m_phase_reference_position[0] -
                                      m_local_phase_reference_position[0]) +
                      direction[1] * (antenna->m_phase_reference_position[1] -
                                      m_local_phase_reference_position[1]) +
                      direction[2] * (antenna->m_phase_reference_position[2] -
                                      m_local_phase_reference_position[2]);

    double phase = -2 * M_PI * dl / (common::c / freq);
    result.push_back({std::sin(phase), std::cos(phase)});
  }
  return result;
}

std::vector<std::pair<std::complex<double>, std::complex<double>>>
BeamFormer::compute_weights(const vector3r_t &pointing, double freq) const {
  // Get geometric response for pointing direction
  auto geometric_response = compute_geometric_response(freq, pointing);

  // Initialize and fill result
  double weight_sum[2] = {0.0, 0.0};
  std::vector<std::pair<std::complex<double>, std::complex<double>>> result(
      geometric_response.size());
  for (std::size_t antenna_idx = 0; antenna_idx < m_antennas.size();
       ++antenna_idx) {
    // Compute conjugate of geometric response
    auto phasor_conj = std::conj(geometric_response[antenna_idx]);
    // Compute the delays in x/y direction
    result.push_back(
        {phasor_conj * (1.0 * m_antennas[antenna_idx]->m_enabled[0]),
         phasor_conj * (1.0 * m_antennas[antenna_idx]->m_enabled[1])});
    weight_sum[0] += (1.0 * m_antennas[antenna_idx]->m_enabled[0]);
    weight_sum[1] += (1.0 * m_antennas[antenna_idx]->m_enabled[1]);
  }

  // Normalize the weight by the number of antennas
  for (std::size_t antenna_idx = 0; antenna_idx < m_antennas.size();
       ++antenna_idx) {
    result[antenna_idx].first /= weight_sum[0];
    result[antenna_idx].second /= weight_sum[1];
  }

  return result;
}

matrix22c_t BeamFormer::local_response(real_t time, real_t freq,
                                       const vector3r_t &direction,
                                       const Options &options) const {
  // Weights based on pointing direction of beam
  auto weights = compute_weights(options.station0, options.freq0);
  // Weights based on direction of interest
  auto geometric_response = compute_geometric_response(freq, direction);

  matrix22c_t result = {0};
  for (unsigned int antenna_idx = 0; antenna_idx < m_antennas.size();
       ++antenna_idx) {
    auto antenna = m_antennas[antenna_idx];
    auto antenna_weight = weights[antenna_idx];
    auto antenna_geometric_reponse = geometric_response[antenna_idx];

    matrix22c_t antenna_response =
        antenna->response(time, freq, direction, options);

    result[0][0] += antenna_weight.first * antenna_geometric_reponse *
                    antenna_response[0][0];
    result[0][1] += antenna_weight.first * antenna_geometric_reponse *
                    antenna_response[0][1];
    result[1][0] += antenna_weight.second * antenna_geometric_reponse *
                    antenna_response[1][0];
    result[1][1] += antenna_weight.second * antenna_geometric_reponse *
                    antenna_response[1][1];
  }
  return result;
}
}  // namespace everybeam
