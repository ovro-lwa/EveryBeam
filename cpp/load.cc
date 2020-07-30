#include "load.h"
#include "options.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <stdexcept>

using namespace everybeam;
namespace everybeam {
namespace {
enum TelescopeType {
  kUnknownTelescope,
  kLofarTelescope,
  kVLATelescope,
  kATCATelescope,
  kMWATelescope
};

TelescopeType Convert(const std::string &telescope_name) {
  if (telescope_name == "LOFAR") return kLofarTelescope;
  // Maybe more elegant with boost::to_upper_copy()?
  else if (telescope_name.compare(0, 4, "EVLA") == 0)
    return kVLATelescope;
  else if (telescope_name.compare(0, 4, "ATCA") == 0)
    return kATCATelescope;
  else if (telescope_name == "MWA")
    return kMWATelescope;
  else
    return kUnknownTelescope;
}
}  // namespace

std::unique_ptr<telescope::Telescope> Load(casacore::MeasurementSet &ms,
                                           const Options &options,
                                           const ElementResponseModel model) {
  // Read Telescope name and convert to enum
  casacore::ScalarColumn<casacore::String> telescope_name_col(ms.observation(),
                                                              "TELESCOPE_NAME");

  TelescopeType telescope_name = Convert(telescope_name_col(0));
  switch (telescope_name) {
    case kLofarTelescope: {
      std::unique_ptr<telescope::Telescope> telescope =
          std::unique_ptr<telescope::Telescope>(
              new telescope::LOFAR(ms, model, options));
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
    default:
      std::stringstream message;
      message << "The requested telescope type " << telescope_name_col(0)
              << " is not implemented.";
      throw std::runtime_error(message.str());
  }
};
}  // namespace everybeam