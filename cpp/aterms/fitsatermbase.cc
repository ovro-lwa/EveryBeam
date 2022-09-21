// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "fitsatermbase.h"
#include "atermresampler.h"

#include <algorithm>
#include <limits>

namespace everybeam {
namespace aterms {

FitsATermBase::FitsATermBase(
    size_t n_antennas, const aocommon::CoordinateSystem& coordinate_system,
    size_t max_support)
    : cache_(n_antennas * 4 * coordinate_system.width *
             coordinate_system.height),
      cur_timeindex_(0),
      cur_frequency_(0),
      n_frequencies_(0),
      n_antennas_(n_antennas),
      coordinate_system_(coordinate_system),
      resampler_(coordinate_system, max_support) {}

FitsATermBase::~FitsATermBase() = default;

void FitsATermBase::InitializeFromFiles(
    std::vector<aocommon::FitsReader>& readers) {
  // Sort the readers on observation time
  std::sort(
      readers.begin(), readers.end(),
      [](const aocommon::FitsReader& a, const aocommon::FitsReader& b) -> bool {
        return a.TimeDimensionStart() < b.TimeDimensionStart();
      });
  n_frequencies_ = readers.empty() ? 0 : readers.front().NFrequencies();
  for (size_t reader_index = 0; reader_index != readers.size();
       ++reader_index) {
    const aocommon::FitsReader& reader = readers[reader_index];
    if (n_frequencies_ != reader.NFrequencies()) {
      throw std::runtime_error(
          "A-term FITS files have inconsistent number of frequencies");
    }
    if (reader.NAntennas() != n_antennas_) {
      std::ostringstream str;
      str << "FITS file for A-terms has incorrect number of antennas. "
             "Measurement set has "
          << n_antennas_ << " antennas, a-term FITS file has "
          << reader.NAntennas() << " antennas.";
      throw std::runtime_error(str.str());
    }
    const double time0 = reader.TimeDimensionStart();
    if (!timesteps_.empty() && time0 < timesteps_.back().time) {
      throw std::runtime_error(
          "Time axis of FITS files seem to overlap (start of fitsfile " +
          std::to_string(reader_index) + " (t=" + std::to_string(time0) +
          " was before end of previous fitsfile)");
    }
    for (size_t i = 0; i != reader.NTimesteps(); ++i) {
      timesteps_.emplace_back(time0 + i * reader.TimeDimensionIncr(),
                              reader_index, i);
    }
  }
  cur_timeindex_ = std::numeric_limits<size_t>::max();
  cur_frequency_ = std::numeric_limits<double>::max();
}

double FitsATermBase::AverageUpdateTime() const {
  if (timesteps_.size() < 2) {
    return 60.0 * 30.0;
  } else {
    return (timesteps_.back().time - timesteps_.front().time) /
           (timesteps_.size() - 1);
  }
}

bool FitsATermBase::FindFilePosition(std::complex<float>* buffer, double time,
                                     double frequency, size_t& timeindex,
                                     bool& requires_recalculation) {
  requires_recalculation = false;
  if (cur_timeindex_ == std::numeric_limits<size_t>::max()) {
    requires_recalculation = true;
    cache_.Reset();
    cur_timeindex_ = 0;
  }

  bool finished_search = false;
  while (cur_timeindex_ + 1 < timesteps_.size() && !finished_search) {
    // Do we need to calculate a next timestep?
    double current_time = timesteps_[cur_timeindex_].time;
    double next_time = timesteps_[cur_timeindex_ + 1].time;
    // If we are closer to the next timestep, use the next.
    if (std::fabs(next_time - time) < std::fabs(current_time - time)) {
      ++cur_timeindex_;
      requires_recalculation = true;
      cache_.Reset();
      finished_search = false;
    } else {
      finished_search = true;
    }
  }
  timeindex = cur_timeindex_;

  if (!requires_recalculation) {
    // If we are here, it means that the timestep didn't
    // change. So if the frequency also didn't change, we're done...
    if (cur_frequency_ == frequency) return false;
    // If it did change: do we have this frequency in the cache?
    size_t cache_index = cache_.Find(frequency);
    if (cache_index == Cache::kNotFound) {
      requires_recalculation = true;
    } else {
      cache_.Get(cache_index, buffer);
      cur_frequency_ = frequency;
      return true;
    }
  }

  return true;
}

void FitsATermBase::StoreInCache(double frequency,
                                 const std::complex<float>* buffer) {
  cur_frequency_ = frequency;
  cache_.Store(frequency, buffer);
}

}  // namespace aterms
}  // namespace everybeam
