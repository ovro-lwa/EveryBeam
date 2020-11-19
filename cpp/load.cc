// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "load.h"

#include "telescope/lofar.h"
#include "telescope/dish.h"
#include "telescope/mwa.h"
#include "telescope/oskar.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <string>
#include <stdexcept>

namespace everybeam {
TelescopeType GetTelescopeType(const casacore::MeasurementSet &ms) {
  // Read Telescope name and convert to enum
  casacore::ScalarColumn<casacore::String> telescope_name_col(ms.observation(),
                                                              "TELESCOPE_NAME");
  std::string telescope_name = telescope_name_col(0);
  std::for_each(telescope_name.begin(), telescope_name.end(),
                [](char &c) { c = ::toupper(c); });

  if (telescope_name == "LOFAR")
    return kLofarTelescope;
  else if (telescope_name == "AARTFAAC")
    return kAARTFAAC;
  else if (telescope_name.compare(0, 4, "EVLA") == 0)
    return kVLATelescope;
  else if (telescope_name.compare(0, 4, "ATCA") == 0)
    return kATCATelescope;
  else if (telescope_name == "MWA")
    return kMWATelescope;
  // check if telescope_name starts with "OSKAR"
  else if (telescope_name.rfind("OSKAR", 0) == 0)
    return kOSKARTelescope;
  else
    return kUnknownTelescope;
}

std::unique_ptr<telescope::Telescope> Load(const casacore::MeasurementSet &ms,
                                           const Options &options) {
  std::unique_ptr<telescope::Telescope> telescope;
  const TelescopeType telescope_name = GetTelescopeType(ms);
  switch (telescope_name) {
    case kAARTFAAC:  // fall through
    case kLofarTelescope:
      telescope.reset(new telescope::LOFAR(ms, options));
      break;
    case kATCATelescope:  // fall through
    case kVLATelescope:
      telescope.reset(new telescope::Dish(ms, options));
      break;
    case kMWATelescope: {
      telescope.reset(new telescope::MWA(ms, options));
      break;
    }
    case kOSKARTelescope: {
      telescope.reset(new telescope::OSKAR(ms, options));
      break;
    }
    default:
      casacore::ScalarColumn<casacore::String> telescope_name_col(
          ms.observation(), "TELESCOPE_NAME");
      std::stringstream message;
      message << "The requested telescope type " << telescope_name_col(0)
              << " is not implemented.";
      throw std::runtime_error(message.str());
  }
  return telescope;
}

std::unique_ptr<telescope::Telescope> Load(const std::string &ms_name,
                                           const Options &options) {
  casacore::MeasurementSet ms(ms_name);
  return Load(ms, options);
}
}  // namespace everybeam
