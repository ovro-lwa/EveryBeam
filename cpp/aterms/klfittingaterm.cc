// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "klfittingaterm.h"

#include <aocommon/imagecoordinates.h>

#include "klfitter.h"

using schaapcommon::h5parm::AxisInfo;
using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::SolTab;

namespace everybeam {
namespace aterms {

KlFittingATerm::KlFittingATerm(
    const std::vector<std::string>& station_names_ms,
    const aocommon::CoordinateSystem& coordinate_system, int order,
    bool use_phasor_fit)
    : station_names_ms_(station_names_ms),
      coordinate_system_(coordinate_system),
      order_(order),
      phase_soltab_(),
      update_interval_(0),
      kl_fitter_(),
      last_aterm_update_(-1),
      nr_directions_(0),
      use_phasor_fit_(use_phasor_fit) {
  if (coordinate_system.height != coordinate_system.width) {
    throw std::runtime_error(
        "Non-square aterms not supported by KlFittingATerm.");
  }
}

// Default destructor defined here, because here the size of a KlFitting
// object is known
KlFittingATerm::~KlFittingATerm() = default;

void KlFittingATerm::Open(const std::string& filename) {
  H5Parm h5parmfile(filename);
  phase_soltab_ = h5parmfile.GetSolTab("phase000");
  const schaapcommon::h5parm::AxisInfo dir_axis = phase_soltab_.GetAxis("dir");
  nr_directions_ = dir_axis.size;
  const std::vector<schaapcommon::h5parm::H5Parm::source_t> sources =
      h5parmfile.GetSources();
  const std::vector<std::string> directions =
      phase_soltab_.GetStringAxis("dir");

  // If no station selection is given, select all stations in the h5parm file
  if (station_names_ms_.empty()) {
    station_names_ms_ = phase_soltab_.GetStringAxis("ant");
  }

  std::map<std::string, int> source_map;
  for (size_t i = 0; i < sources.size(); ++i) {
    source_map[sources[i].name] = i;
  }

  // Source directions in floating point pixel coordinates
  // Order of the coordinates matches the order of the axes in the output array
  // The fastest changing index in the output array (last one in C(++) / row
  // major ordering) corresponds to the ra/l/x axis, so x comes last here.
  std::vector<std::pair<float, float>> directions_yx;

  for (const std::string& dir : directions) {
    const int source_idx = source_map[dir];
    const float ra = sources[source_idx].dir[0];
    const float dec = sources[source_idx].dir[1];
    float l, m;
    aocommon::ImageCoordinates::RaDecToLM(ra, dec, float(coordinate_system_.ra),
                                          float(coordinate_system_.dec), l, m);
    float x, y;
    aocommon::ImageCoordinates::LMToXYfloat(
        float(l - coordinate_system_.l_shift),
        float(m - coordinate_system_.m_shift),
        -float(coordinate_system_
                   .dl),  // LMToXYfloat already applies a minus sign to the
                          // l-axis. For the usual l,m coordinate system, dl is
                          // negative These two minuses cancel eachother,
                          // therefore another minus is needed here
        float(coordinate_system_.dm), coordinate_system_.width,
        coordinate_system_.height, x, y);
    directions_yx.emplace_back(y, x);
  }

  const int subgrid_size = coordinate_system_.height;

  kl_fitter_ = std::make_unique<KlFitter>(subgrid_size, order_, directions_yx);
}

bool KlFittingATerm::Calculate(std::complex<float>* buffer, double time,
                               double frequency,
                               [[maybe_unused]] size_t field_id,
                               [[maybe_unused]] const double* uvw_in_m) {
  const bool outdated = std::fabs(time - last_aterm_update_) > update_interval_;
  if (!outdated) return false;
  last_aterm_update_ = time;

  const size_t time_index = phase_soltab_.GetTimeIndex(time);
  const size_t freq_index = phase_soltab_.GetFreqIndex(frequency);

  std::vector<float> solutions_ref(nr_directions_, 0.0);
  int ant = 0;
  for (auto const& station_name : station_names_ms_) {
    const size_t n_times = 1;
    const size_t timestep = 1;

    const size_t n_freq = 1, freq_step = 1;
    const size_t pol = 0;

    std::vector<float> solutions;
    solutions.reserve(nr_directions_);

    for (size_t dir = 0; dir < nr_directions_; ++dir) {
      double phase =
          phase_soltab_.GetValues(station_name, time_index, n_times, timestep,
                                  freq_index, n_freq, freq_step, pol, dir)[0];
      solutions.push_back(phase - solutions_ref[dir]);
    }

    if (ant == 0) {
      std::swap(solutions, solutions_ref);
    }

    const int subgrid_size = coordinate_system_.height;

    if (use_phasor_fit_) {
      std::vector<float> screen_real(subgrid_size * subgrid_size);
      std::vector<float> screen_imag(subgrid_size * subgrid_size);
      std::vector<float> solutions_phasor_real;
      std::vector<float> solutions_phasor_imag;
      std::transform(solutions.begin(), solutions.end(),
                     back_inserter(solutions_phasor_real), cos);
      std::transform(solutions.begin(), solutions.end(),
                     back_inserter(solutions_phasor_imag), sin);

      kl_fitter_->Evaluate(solutions_phasor_real, screen_real.data());
      kl_fitter_->Evaluate(solutions_phasor_imag, screen_imag.data());

      for (size_t i = 0; i < (subgrid_size * subgrid_size); ++i) {
        // Store output on "diagonal" of Jones matrix
        std::complex<float> phasor(
            std::complex<double>(screen_real[i], screen_imag[i]));
        phasor /= abs(phasor);
        buffer[0] = phasor;
        buffer[1] = 0.0;
        buffer[2] = 0.0;
        buffer[3] = phasor;
        buffer += 4;
      }
    } else {
      std::vector<float> phase_screen(subgrid_size * subgrid_size);
      kl_fitter_->Evaluate(solutions, phase_screen.data());
      for (size_t i = 0; i < (subgrid_size * subgrid_size); ++i) {
        // Store output on "diagonal" of Jones matrix
        buffer[0] =
            std::complex<float>(exp(std::complex<double>(0, phase_screen[i])));
        buffer[1] = 0.0;
        buffer[2] = 0.0;
        buffer[3] =
            std::complex<float>(exp(std::complex<double>(0, phase_screen[i])));
        buffer += 4;
      }
    }
    ant++;
  }
  return true;
}

}  // namespace aterms
}  // namespace everybeam
