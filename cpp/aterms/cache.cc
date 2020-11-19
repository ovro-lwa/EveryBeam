// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "cache.h"

#include <limits>

using everybeam::aterms::Cache;
const size_t Cache::kNotFound = std::numeric_limits<size_t>::max();