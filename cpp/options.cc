// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "options.h"

#include "config.h"

#include <cstdlib>
#include <cstring>

namespace everybeam {

// TODO(RAP-260) Improve this path lookup.
std::filesystem::path GetDataDirectory() {
  const char* envvar;
  if (std::strcmp(EVERYBEAM_DATADIR, EVERYBEAM_FULL_DATADIR) == 0)
    return std::filesystem::path(EVERYBEAM_DATADIR);
  if ((envvar = std::getenv("EVERYBEAM_DATADIR")))
    return std::filesystem::path(envvar);
  if ((envvar = std::getenv("CONDA_PREFIX")))
    return std::filesystem::path(envvar) / EVERYBEAM_DATADIR;
  if ((envvar = std::getenv("VIRTUAL_ENV")))
    return std::filesystem::path(envvar) / EVERYBEAM_DATADIR;
  return std::filesystem::path(EVERYBEAM_FULL_DATADIR);
}
}  // namespace everybeam
