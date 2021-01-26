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

    ampl_polynomial_ = std::unique_ptr<LagrangePolynomial>(
        new LagrangePolynomial(amplitude_soltab_.back().GetAxis("dir").size));
    phase_polynomial_ = std::unique_ptr<LagrangePolynomial>(
        new LagrangePolynomial(phase_soltab_.back().GetAxis("dir").size));

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
      if (station_names_ms_[i] != station_names_ampl[i])
        throw std::runtime_error(
            "At index " + std::to_string(i) + ": station name " +
            station_names_ms_[i] +
            " provided by ms does not match station name " +
            station_names_ampl[i] + " provided by amplitude soltab");

      if (station_names_ms_[i] != station_names_phase[i])
        throw std::runtime_error(
            "At index " + std::to_string(i) + ": station name " +
            station_names_ms_[i] +
            " provided by ms does not match station name " +
            station_names_phase[i] + " provided by phase soltab");
    }
  }
}

bool H5ParmATerm::Calculate(std::complex<float>* buffer, double time,
                            double frequency, size_t, const double*) {
  bool outdated = std::fabs(time - last_aterm_update_) > update_interval_;
  if (!outdated) return false;
  last_aterm_update_ = time;

  hsize_t tindex_ampl = amplitude_soltab_[0].GetTimeIndex(time);
  hsize_t tindex_phase = phase_soltab_[0].GetTimeIndex(time);
  bool recalculate_amplitude = (tindex_ampl != last_ampl_index_);
  bool recalculate_phase = (tindex_phase != last_phase_index_);

  // Outer loop may be over the y coordinates when implementing multi-threading
  // for (const auto& name : station_names_ms_) {
  for (size_t i = 0; i < station_names_ms_.size(); ++i) {
    std::string name = station_names_ms_[i];
    size_t station_offset =
        i * coordinate_system_.height * coordinate_system_.width;
    for (size_t y = 0; y < coordinate_system_.height; ++y) {
      for (size_t x = 0; x < coordinate_system_.width; ++x) {
        double l, m;
        aocommon::ImageCoordinates::XYToLM(
            x, y, coordinate_system_.dl, coordinate_system_.dm,
            coordinate_system_.width, coordinate_system_.height, l, m);

        l += coordinate_system_.phase_centre_dl;
        m += coordinate_system_.phase_centre_dm;

        size_t offset = station_offset + y * coordinate_system_.width + x;
        std::complex<float> output =
            ExpandComplexExp(name, tindex_ampl, tindex_phase, frequency, l, m,
                             recalculate_amplitude, recalculate_phase, offset);

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
  last_ampl_index_ = tindex_ampl;
  last_phase_index_ = tindex_phase;
  return true;
}

void H5ParmATerm::ReadCoeffs(SolTab& soltab, const std::string& station_name,
                             std::vector<float>& coeffs, hsize_t tindex,
                             double) {
  size_t ntime = 1;
  size_t timestep = 1;

  // Not yet relevant
  size_t freq_start = 0, nfreq = 1, freqstep = 1;
  size_t pol = 0;

  for (size_t idx = 0; idx < coeffs.size(); ++idx) {
    coeffs[idx] = soltab.GetValues(station_name, tindex, ntime, timestep,
                                   freq_start, nfreq, freqstep, pol, idx)[0];
  }
}

std::complex<float> H5ParmATerm::ExpandComplexExp(
    const std::string& station_name, hsize_t ampl_tindex, hsize_t phase_tindex,
    double freq, double l, double m, bool recalculate_ampl,
    bool recalculate_phase, size_t offset) {
  // Expand amplitude coeffs
  if (recalculate_ampl) {
    std::vector<float> ampl_coeffs(ampl_polynomial_->GetNrCoeffs());
    ReadCoeffs(amplitude_soltab_[0], station_name, ampl_coeffs, ampl_tindex,
               freq);
    amplitude_cache_[offset] = ampl_polynomial_->Evaluate(l, m, ampl_coeffs);
  }

  // Expand phase coeffs
  if (recalculate_phase) {
    std::vector<float> phase_coeffs(phase_polynomial_->GetNrCoeffs());
    ReadCoeffs(phase_soltab_[0], station_name, phase_coeffs, phase_tindex,
               freq);
    phase_cache_[offset] = phase_polynomial_->Evaluate(l, m, phase_coeffs);
  }
  return std::complex<float>{amplitude_cache_[offset] *
                             exp(std::complex<float>(0, phase_cache_[offset]))};
}
}  // namespace aterms
}  // namespace everybeam