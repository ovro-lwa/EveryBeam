#include "fourierfittingaterm.h"

#include <cmath>
#include <algorithm>
#include <iostream>

#include <aocommon/imagecoordinates.h>

#include "fourierfitter.h"

using schaapcommon::h5parm::AxisInfo;
using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::SolTab;

namespace everybeam {
namespace aterms {

FourierFittingATerm::FourierFittingATerm(
    const std::vector<std::string>& station_names_ms,
    const aocommon::CoordinateSystem& coordinate_system, int support)
    : fourier_fitter_(),
      station_names_ms_(station_names_ms),
      coordinate_system_(coordinate_system),
      support_(support),
      update_interval_(0),
      last_aterm_update_(-1) {
  if (coordinate_system.height != coordinate_system.width) {
    throw std::runtime_error(
        "Non-square aterms not supported by FourierFittingATerm.");
  }
}

// Default destructor defined here, because here the size of a FourierFitting
// object is known
FourierFittingATerm::~FourierFittingATerm() = default;

void FourierFittingATerm::Open(const std::string& filename) {
  H5Parm h5parmfile(filename);
  phase_soltab_ = h5parmfile.GetSolTab("phase000");
  const schaapcommon::h5parm::AxisInfo dir_axis = phase_soltab_.GetAxis("dir");
  nr_directions_ = dir_axis.size;
  const std::vector<schaapcommon::h5parm::H5Parm::source_t> sources =
      h5parmfile.GetSources();
  const std::vector<std::string> directions =
      phase_soltab_.GetStringAxis("dir");
  std::map<std::string, int> source_map;
  for (size_t i = 0; i < sources.size(); ++i) {
    source_map[sources[i].name] = i;
  }
  std::vector<std::pair<float, float>> directions_xy;
  for (const std::string& dir : directions) {
    const int source_idx = source_map[dir];
    const float ra = sources[source_idx].dir[0];
    const float dec = sources[source_idx].dir[1];
    float l, m;
    aocommon::ImageCoordinates::RaDecToLM(ra, dec, float(coordinate_system_.ra),
                                          float(coordinate_system_.dec), l, m);
    float x, y;
    aocommon::ImageCoordinates::LMToXYfloat(
        l, m, float(coordinate_system_.dl), float(coordinate_system_.dm),
        coordinate_system_.width, coordinate_system_.height, x, y);
    directions_xy.emplace_back(x, y);
  }

  const int subgrid_size = coordinate_system_.height;

  fourier_fitter_ =
      std::make_unique<FourierFitter>(subgrid_size, support_, directions_xy);
}

bool FourierFittingATerm::Calculate(std::complex<float>* buffer, double time,
                                    [[maybe_unused]] double frequency,
                                    [[maybe_unused]] size_t field_id,
                                    [[maybe_unused]] const double* uvm_in_m) {
  const bool outdated = std::fabs(time - last_aterm_update_) > update_interval_;
  if (!outdated) return false;
  last_aterm_update_ = time;

  const size_t time_index = phase_soltab_.GetTimeIndex(time);
  const size_t freq_index = phase_soltab_.GetFreqIndex(frequency);

  for (auto const& station_name : station_names_ms_) {
    const size_t n_times = 1;
    const size_t timestep = 1;

    const size_t n_freq = 1, freq_step = 1;
    const size_t pol = 0;

    std::vector<std::complex<float>> solutions;
    solutions.reserve(nr_directions_);

    for (size_t dir = 0; dir < nr_directions_; ++dir) {
      double phase =
          phase_soltab_.GetValues(station_name, time_index, n_times, timestep,
                                  freq_index, n_freq, freq_step, pol, dir)[0];
      solutions.push_back(
          std::complex<float>(exp(std::complex<double>(0, phase))));
    }
    const int subgrid_size = coordinate_system_.height;
    std::vector<std::complex<float>> screen(subgrid_size * subgrid_size);
    fourier_fitter_->Evaluate(solutions, screen.data());

    for (size_t i = 0; i < (subgrid_size * subgrid_size); ++i) {
      // Store output on "diagonal" of Jones matrix
      buffer[0] = screen[i];
      buffer[1] = 0.0;
      buffer[2] = 0.0;
      buffer[3] = screen[i];
      buffer += 4;
    }
  }
  return true;
}

}  // namespace aterms
}  // namespace everybeam
