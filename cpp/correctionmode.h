// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_CORRECTIONMODE_H_
#define EVERYBEAM_CORRECTIONMODE_H_

#include <boost/algorithm/string/case_conv.hpp>

#include <string>

namespace everybeam {
/**
 * Describes to what level the beam is corrected: array, element or both.
 * These may have different meanings for different telescopes. For LOFAR HBA,
 * the element beam is the dipole beam, and the array beam is the combination
 * of the tile factor and station-array factor.
 */
enum class CorrectionMode { kNone, kFull, kArrayFactor, kElement };

inline std::string ToString(CorrectionMode mode) {
  switch (mode) {
    case CorrectionMode::kNone:
      return "None";
    case CorrectionMode::kFull:
      return "Full";
    case CorrectionMode::kArrayFactor:
      return "ArrayFactor";
    case CorrectionMode::kElement:
      return "Element";
  }
  throw std::runtime_error("Invalid correction mode");
}

inline CorrectionMode ParseCorrectionMode(const std::string& str) {
  const std::string lower_str = boost::algorithm::to_lower_copy(str);
  if (lower_str == "none")
    return CorrectionMode::kNone;
  else if (lower_str == "full" || lower_str == "default")
    return CorrectionMode::kFull;
  else if (lower_str == "arrayfactor" || lower_str == "array_factor")
    return CorrectionMode::kArrayFactor;
  else if (lower_str == "element")
    return CorrectionMode::kElement;
  else
    throw std::runtime_error(
        "Invalid beam correction mode \'" + str +
        "\', options are: None, Default, Full, ArrayFactor or Element");
}
}  // namespace everybeam

#endif
