// Options.h: Class for specifying the telescope response options.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_OPTIONS_H_
#define EVERYBEAM_OPTIONS_H_

#include <filesystem>
#include <string>
#include <vector>

#include "elementresponse.h"
#include "beammode.h"
#include "beamnormalisationmode.h"

namespace everybeam {

/**
 * @brief Class/Struct specifying everybeam Options. Needs further
 * implementation!
 *
 */
struct Options {
  // Path to (MWA/LOBES) coefficients file:
  // - in case of MWA, this should be the path to the file (absolute/relative)
  // - in case of LOBES, this should be the directory containing LOBES
  // coefficient files
  std::string coeff_path = "";
  BeamNormalisationMode beam_normalisation_mode = BeamNormalisationMode::kNone;

  // LOFAR specific
  bool use_channel_frequency = true;
  std::string data_column_name = "DATA";
  ElementResponseModel element_response_model = ElementResponseModel::kHamaker;
  BeamMode beam_mode = BeamMode::kFull;

  // MWA specific (Lofar probably will follow)
  bool frequency_interpolation = false;
};

struct ATermSettings {
  // Path to coefficients file
  std::string coeff_path = ".";
  // Save aterm fits files?
  bool save_aterms = false;
  // Prefix for the aterm fits files
  std::string save_aterms_prefix = "wsclean";
  // Default for the data column name
  std::string data_column_name = "DATA";

  std::vector<std::string> filenames = std::vector<std::string>();
  size_t max_support = 32;

  // Time interval (s) before a new aterm will be computed
  double aterm_update_interval = 300.;
  size_t padded_image_width = 0, padded_image_height = 0,
         trimmed_image_width = 0, trimmed_image_height = 0;
};

/**
 * Get the path prefix to the EveryBeam data directory.
 * Which path is returned depends on your environment, and the way that
 * @c EveryBeam was built. If at CMake configure-time an absolute path was
 * given to @c EVERYBEAM_DATADIR, then this absolute path will be returned
 * as prefix. Otherwise, the following @a environment variables will be
 * checked, in the given order:
 * - @c EVERYBEAM_DATADIR: path prefix set by a user;
 * - @c CONDA_PREFIX: set when running in a Conda environment;
 * - @c VIRTUAL_ENV: set when running inside a Python virtual environment
 * If none of these environment variables are set, then the full path that
 * was set at compile-time is returned.
 * @return Path prefix to the EveryBeam data directory.
 */
std::filesystem::path GetDataDirectory();

}  // namespace everybeam
#endif  // EVERYBEAM_OPTIONS_H_
