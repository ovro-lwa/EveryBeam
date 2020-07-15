#include "load.h"
#include "options.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>

using namespace everybeam;
namespace everybeam {
namespace {
enum TelescopeType { kUnknownTelescope, kLofarTelescope };

TelescopeType Convert(const std::string &str) {
  if (str == "LOFAR")
    return kLofarTelescope;
  else
    return kUnknownTelescope;
}
}  // namespace

std::unique_ptr<telescope::Telescope> Load(casacore::MeasurementSet &ms,
                                           const ElementResponseModel model,
                                           const Options &options) {
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
    default:
      std::stringstream message;
      message << "The requested telescope type " << telescope_name_col(0)
              << " is not implemented.";
      throw std::runtime_error(message.str());
  }
};
}  // namespace everybeam