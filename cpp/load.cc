#include "load.h"

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
  else if (str.rfind("OSKAR", 0) == 0)
    return kOskarTelescope;
  else
    return kUnknownTelescope;
}

std::unique_ptr<telescope::Telescope> Load(casacore::MeasurementSet &ms,
                                           const Options &options) {
  TelescopeType telescope_name = GetTelescopeType(ms);
  switch (telescope_name) {
    case kAARTFAAC:
    case kLofarTelescope: {
      std::unique_ptr<telescope::Telescope> telescope =
          std::unique_ptr<telescope::Telescope>(
              new telescope::LOFAR(ms, options));
      return telescope;
    }
    case kATCATelescope: {
      std::unique_ptr<telescope::Telescope> telescope =
          std::unique_ptr<telescope::Telescope>(
              new telescope::Dish(ms, options));
      return telescope;
    }
    case kVLATelescope: {
      std::unique_ptr<telescope::Telescope> telescope =
          std::unique_ptr<telescope::Telescope>(
              new telescope::Dish(ms, options));
      return telescope;
    }
    case kMWATelescope: {
      std::unique_ptr<telescope::Telescope> telescope =
          std::unique_ptr<telescope::Telescope>(
              new telescope::MWA(ms, options));
      return telescope;
    }
    case kOskarTelescope: {
      std::unique_ptr<telescope::Telescope> telescope =
          std::unique_ptr<telescope::Telescope>(
              new telescope::OSKAR(ms, model, options));
      return telescope;
    }
    default:
      casacore::ScalarColumn<casacore::String> telescope_name_col(
          ms.observation(), "TELESCOPE_NAME");
      std::stringstream message;
      message << "The requested telescope type " << telescope_name_col(0)
              << " is not implemented.";
      throw std::runtime_error(message.str());
  }
};
}  // namespace everybeam
