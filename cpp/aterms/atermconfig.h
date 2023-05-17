// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_ATERM_CONFIG_H
#define EVERYBEAM_ATERMS_ATERM_CONFIG_H

#include "atermbase.h"
#include "atermbeam.h"
#include "../options.h"

#include <string>
#include <vector>

#include <aocommon/coordinatesystem.h>
#include <aocommon/uvector.h>

namespace casacore {
class MeasurementSet;
}

namespace everybeam {
namespace aterms {

class ParsetProvider;

class [[gnu::visibility("default")]] ATermConfig final : public ATermBase {
 public:
  ATermConfig(size_t n_antennas,
              const aocommon::CoordinateSystem& coordinate_system,
              const everybeam::ATermSettings& settings)
      : n_antennas_(n_antennas),
        coordinate_system_(coordinate_system),
        settings_(settings),
        aterms_(),
        previous_aterm_values_() {}

  /**
   * @brief Read parset information
   *
   * @param ms Measurement set
   * @param reader Implementation of ParsetProvider
   * @param ms_filename ms filename, only relevant for paf-type aterm
   */
  void Read(const casacore::MeasurementSet& ms, const ParsetProvider& reader,
            const std::string& ms_filename = "");

  /** Reimplemented from ATermBase */
  bool Calculate(std::complex<float> * buffer, double time, double frequency,
                 size_t fieldId, const double* uvwInM) override;

  /** Reimplemented from ATermBase */
  double AverageUpdateTime() const override {
    double avgTime = aterms_.front()->AverageUpdateTime();
    for (size_t i = 1; i < aterms_.size(); ++i)
      avgTime = std::min(avgTime, aterms_[i]->AverageUpdateTime());
    return avgTime;
  }

  /**
   * @brief Static method for constructing an (EveryBeam)ATerm
   *
   * @param ms Measurement set
   * @param coordinate_system struct with image settings
   * @param settings aterm specific settings
   * @param frequency_interpolation Interpolate between frequencies? Relevant
   * for MWA only.
   * @param beam_normalisation_mode What type of normalisation will be applied?
   * @param use_channel_frequency Use channel frequency
   * @param element_response_model Element response model
   * @return std::unique_ptr<ATermBeam>
   */
  static std::unique_ptr<ATermBeam> GetATermBeam(
      const casacore::MeasurementSet& ms,
      const aocommon::CoordinateSystem& coordinate_system,
      const ATermSettings& settings, bool frequency_interpolation,
      const std::string& beam_normalisation_mode, bool use_channel_frequency,
      const std::string& element_response_model,
      const std::string& beam_mode = "full");

  /**
   * @brief Static method to construct an everybeam::Options struct
   * from user-settings
   *
   * @param ms Measurement set
   * @param settings ATermSettings
   * @param frequency_interpolation Interpolate between frequencies? Relevant
   * for MWA only.
   * @param beam_normalisation_mode What type of normalisation will be applied?
   * @param use_channel_frequency Use channel frequency?
   * @param element_response_model Element response model
   * @return everybeam::Options
   */
  static everybeam::Options ConvertToEBOptions(
      const casacore::MeasurementSet& ms, const ATermSettings& settings,
      bool frequency_interpolation, const std::string& beam_normalisation_mode,
      bool use_channel_frequency, const std::string& element_response_model,
      const std::string& beam_mode = "full");

 private:
  const size_t n_antennas_;
  const aocommon::CoordinateSystem coordinate_system_;
  const ATermSettings settings_;
  std::vector<std::unique_ptr<ATermBase>> aterms_;
  std::vector<aocommon::UVector<std::complex<float>>> previous_aterm_values_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
