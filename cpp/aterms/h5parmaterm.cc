#include "h5parmaterm.h"

#include "../common/fftresampler.h"

#include <aocommon/imagecoordinates.h>
#include <cmath>
#include <algorithm>

using schaapcommon::h5parm::AxisInfo;
using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::SolTab;

namespace everybeam {
namespace aterms {

H5ParmATerm::H5ParmATerm(const std::vector<std::string>& station_names_ms,
                         const coords::CoordinateSystem& coordinate_system)
    : station_names_ms_(station_names_ms),
      coordinate_system_(coordinate_system),
      update_interval_(0),
      last_aterm_update_(-1),
      last_ampl_index_(std::numeric_limits<hsize_t>::max()),
      last_phase_index_(std::numeric_limits<hsize_t>::max()),
      amplitude_cache_(station_names_ms_.size() * coordinate_system_.width *
                       coordinate_system_.height),
      phase_cache_(station_names_ms_.size() * coordinate_system_.width *
                   coordinate_system_.height) {}

void H5ParmATerm::Open(const std::vector<std::string>& filenames) {
  if (filenames.size() > 1) {
    throw std::runtime_error("Multiple h5parm input files not (yet) supported");
  }
  for (auto const& filename : filenames) {
    H5Parm h5parmfile(filename);
    // Fill solution tables
    amplitude_soltab_.push_back(h5parmfile.GetSolTab("amplitude_coefficients"));
    phase_soltab_.push_back(h5parmfile.GetSolTab("phase_coefficients"));

    ampl_polynomial_ = std::make_unique<LagrangePolynomial>(
        amplitude_soltab_.back().GetAxis("dir").size);
    phase_polynomial_ = std::make_unique<LagrangePolynomial>(
        phase_soltab_.back().GetAxis("dir").size);

    // Check that antenna names in MS and h5parm match exactly
    std::vector<std::string> station_names_ampl =
        amplitude_soltab_.back().GetStringAxis("ant");
    std::vector<std::string> station_names_phase =
        phase_soltab_.back().GetStringAxis("ant");

    if (station_names_phase.size() != station_names_ampl.size()) {
      throw std::runtime_error(
          "Number of stations in amplitude soltab not equal to number of "
          "stations in phase soltab.");
    }

    if (station_names_ms_.size() != station_names_ampl.size()) {
      throw std::runtime_error(
          "Number of stations in MS not equal to number of stations in h5parm "
          "amplitude soltab.");
    }

    for (size_t i = 0; i < station_names_ms_.size(); ++i) {
      if (station_names_ms_[i] != station_names_ampl[i]) {
        throw std::runtime_error(
            "At index " + std::to_string(i) + ": station name " +
            station_names_ms_[i] +
            " provided by ms does not match station name " +
            station_names_ampl[i] + " provided by amplitude soltab");
      }

      if (station_names_ms_[i] != station_names_phase[i]) {
        throw std::runtime_error(
            "At index " + std::to_string(i) + ": station name " +
            station_names_ms_[i] +
            " provided by ms does not match station name " +
            station_names_phase[i] + " provided by phase soltab");
      }
    }
  }
}

bool H5ParmATerm::Calculate(std::complex<float>* buffer, double time,
                            [[maybe_unused]] double frequency,
                            [[maybe_unused]] size_t field_id,
                            [[maybe_unused]] const double* uvw_in_m) {
  const bool outdated = std::fabs(time - last_aterm_update_) > update_interval_;
  if (!outdated) return false;
  last_aterm_update_ = time;

  const hsize_t time_index_amplitude = amplitude_soltab_[0].GetTimeIndex(time);
  const hsize_t time_index_phase = phase_soltab_[0].GetTimeIndex(time);
  const bool recalculate_amplitude = (time_index_amplitude != last_ampl_index_);
  const bool recalculate_phase = (time_index_phase != last_phase_index_);

  // Initialize placeholders for the y coefficient expansions only once for
  // efficiency reasons
  std::vector<float> scratch_amplitude_coeffs(ampl_polynomial_->GetOrder() + 1);
  std::vector<float> scratch_phase_coeffs(phase_polynomial_->GetOrder() + 1);

  // Outer loop may be over the y coordinates when implementing multi-threading
  for (size_t i = 0; i < station_names_ms_.size(); ++i) {
    const std::string& name = station_names_ms_[i];
    const size_t station_offset =
        i * coordinate_system_.height * coordinate_system_.width;
    for (size_t y = 0; y < coordinate_system_.height; ++y) {
      for (size_t x = 0; x < coordinate_system_.width; ++x) {
        double l, m;
        aocommon::ImageCoordinates::XYToLM(
            x, y, coordinate_system_.dl, coordinate_system_.dm,
            coordinate_system_.width, coordinate_system_.height, l, m);

        l += coordinate_system_.phase_centre_dl;
        m += coordinate_system_.phase_centre_dm;

        const size_t offset = station_offset + y * coordinate_system_.width + x;
        const std::complex<float> output =
            ExpandComplexExp(name, time_index_amplitude, time_index_phase, l, m,
                             recalculate_amplitude, recalculate_phase, offset,
                             scratch_amplitude_coeffs, scratch_phase_coeffs);

        // Store output on "diagonal" of Jones matrix
        buffer[0] = output;
        buffer[1] = 0.0;
        buffer[2] = 0.0;
        buffer[3] = output;
        buffer += 4;
      }
    }
  }
  // Update amplitude and phase index to latest
  last_ampl_index_ = time_index_amplitude;
  last_phase_index_ = time_index_phase;
  return true;
}

void H5ParmATerm::ReadCoeffs(SolTab& soltab, const std::string& station_name,
                             std::vector<float>& coeffs, hsize_t time_index) {
  const size_t n_times = 1;
  const size_t timestep = 1;

  // Not yet relevant
  const size_t freq_start = 0, n_freq = 1, freq_step = 1;
  const size_t pol = 0;

  for (size_t idx = 0; idx < coeffs.size(); ++idx) {
    coeffs[idx] = soltab.GetValues(station_name, time_index, n_times, timestep,
                                   freq_start, n_freq, freq_step, pol, idx)[0];
  }
}

std::complex<float> H5ParmATerm::ExpandComplexExp(
    const std::string& station_name, hsize_t ampl_tindex, hsize_t phase_tindex,
    double l, double m, bool recalculate_ampl, bool recalculate_phase,
    size_t offset, std::vector<float>& scratch_amplitude_coeffs,
    std::vector<float>& scratch_phase_coeffs) {
  if (recalculate_ampl) {
    std::vector<float> ampl_coeffs(ampl_polynomial_->GetNrCoeffs());
    ReadCoeffs(amplitude_soltab_[0], station_name, ampl_coeffs, ampl_tindex);
    amplitude_cache_[offset] =
        ampl_polynomial_->Evaluate(l, m, ampl_coeffs, scratch_amplitude_coeffs);
  }

  if (recalculate_phase) {
    std::vector<float> phase_coeffs(phase_polynomial_->GetNrCoeffs());
    ReadCoeffs(phase_soltab_[0], station_name, phase_coeffs, phase_tindex);
    phase_cache_[offset] =
        phase_polynomial_->Evaluate(l, m, phase_coeffs, scratch_phase_coeffs);
  }
  // Compute complex exponential as Ampl * e^(i*Phase)
  return std::complex<float>{
      amplitude_cache_[offset] *
      std::exp(std::complex<float>(0, phase_cache_[offset]))};
}
}  // namespace aterms
}  // namespace everybeam