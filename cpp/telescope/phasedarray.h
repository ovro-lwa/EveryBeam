// phasedarray.h: Base class for computing the response for phased array
// telescopes (e.g. LOFAR, SKA-LOW)
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_TELESCOPE_PHASEDARRAY_H_
#define EVERYBEAM_TELESCOPE_PHASEDARRAY_H_

#include "../station.h"
#include "../elementresponse.h"
#include "../msreadutils.h"
#include "telescope.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <memory>

namespace everybeam {
namespace telescope {

//! PhasedArray telescope class, is parent to OSKAR and LOFAR
class PhasedArray : public Telescope {
 public:
  /**
   * @brief Construct a new PhasedArray object
   *
   * @param ms MeasurementSet
   * @param model Element Response model
   * @param options telescope options
   */
  PhasedArray(const casacore::MeasurementSet &ms, const Options &options)
      : Telescope(ms, options) {
    stations_.resize(nstations_);
  };

  /**
   * @brief Get station by index
   *
   * @param station_id Station index to retrieve
   * @return std::shared_ptr<Station>
   */
  std::shared_ptr<Station> GetStation(std::size_t station_idx) const {
    // Assert only in DEBUG mode
    assert(station_idx < nstations_);
    return stations_[station_idx];
  }

  // Convenience getters (used in pybindings only)

  //! Get the delay direction
  casacore::MDirection GetDelayDirection() const {
    return ms_properties_.delay_dir;
  }

  //! Get the tile beam direction
  virtual casacore::MDirection GetTileBeamDirection() const {
    return ms_properties_.tile_beam_dir;
  };

  //! Get the preapplied beam direction
  virtual casacore::MDirection GetPreappliedBeamDirection() const {
    return ms_properties_.preapplied_beam_dir;
  };

  //! Get the subband frequency
  double GetSubbandFrequency() const { return ms_properties_.subband_freq; };

  //! Get the number of channels
  size_t GetNrChannels() const { return ms_properties_.channel_count; };

  //! Get the channel frequency given a (zero-based) index
  double GetChannelFrequency(size_t idx) const {
    assert(idx < ms_properties_.channel_count);
    return ms_properties_.channel_freqs[idx];
  };

 protected:
  std::vector<std::shared_ptr<Station>> stations_;

  /**
   * @brief Read stations into vector
   *
   * @param out_it std::vector of stations
   * @param ms measurement set
   * @param model model
   */
  void ReadAllStations(const casacore::MeasurementSet &ms,
                       const ElementResponseModel model) {
    casacore::ROMSAntennaColumns antenna(ms.antenna());

    for (std::size_t i = 0; i < antenna.nrow(); ++i) {
      stations_[i] = ReadStation(ms, i, model);
    }
  };

  std::shared_ptr<Station> ReadStation(const casacore::MeasurementSet &ms,
                                       std::size_t id,
                                       const ElementResponseModel model) const {
    std::shared_ptr<Station> station = ReadSingleStation(ms, id, model);
    return station;
  }

  struct MSProperties {
    double subband_freq;
    casacore::MDirection delay_dir, tile_beam_dir, preapplied_beam_dir;
    size_t channel_count;
    std::vector<double> channel_freqs;
  };

  MSProperties ms_properties_;
};
}  // namespace telescope
}  // namespace everybeam
#endif  // EVERYBEAM_TELESCOPE_PHASEDARRAY_H_
