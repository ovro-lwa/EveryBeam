// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_BEAMNORMALISATIONMODE_H_
#define EVERYBEAM_BEAMNORMALISATIONMODE_H_

#include <boost/algorithm/string/case_conv.hpp>

#include <string>

namespace everybeam {
/**
 * Describes how the beam is normalized:
 *   kNone: unnormalized,
 *   kPreApplied: apply the same correction as has been applied to the
 * visibilities. kAmplitude: scalar amplitude (as in OKSAR) or kFull: full
 * Jones.
 */
enum class BeamNormalisationMode {
  kNone,
  kPreApplied,
  kPreAppliedOrFull,
  kAmplitude,
  kFull
};

inline std::string ToString(BeamNormalisationMode mode) {
  switch (mode) {
    case BeamNormalisationMode::kNone:
      return "None";
    case BeamNormalisationMode::kPreApplied:
      return "PreApplied";
    case BeamNormalisationMode::kPreAppliedOrFull:
      return "PreAppliedOrFull";
    case BeamNormalisationMode::kAmplitude:
      return "Amplitude";
    case BeamNormalisationMode::kFull:
      return "Full";
  }
  throw std::runtime_error("Invalid beam normalisation mode");
}

inline BeamNormalisationMode ParseBeamNormalisationMode(
    const std::string& str) {
  const std::string lower_str = boost::algorithm::to_lower_copy(str);
  if (lower_str == "none")
    return BeamNormalisationMode::kNone;
  else if (lower_str == "preapplied" || lower_str == "pre_applied")
    return BeamNormalisationMode::kPreApplied;
  else if (lower_str == "preappliedorfull" ||
           lower_str == "preapplied_or_full" ||
           lower_str == "pre_applied_or_full")
    return BeamNormalisationMode::kPreAppliedOrFull;
  else if (lower_str == "amplitude")
    return BeamNormalisationMode::kAmplitude;
  else if (lower_str == "full")
    return BeamNormalisationMode::kFull;
  else
    throw std::runtime_error(
        "Invalid beam normalisation mode \'" + str +
        "\', options are: None, PreApplied, Amplitude or Full");
}
}  // namespace everybeam

#endif
