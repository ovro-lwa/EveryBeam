#include "load.h"
#include "options.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>

using namespace everybeam;
namespace everybeam {
namespace {
enum TelescopeType { UNKNOWN_TELESCOPE, LOFAR_TELESCOPE };

TelescopeType Convert(const std::string &str) {
  if (str == "LOFAR")
    return LOFAR_TELESCOPE;
  else
    return UNKNOWN_TELESCOPE;
}

// Default options
// Options defaultOptions;
}  // namespace

telescope::Telescope::Ptr Load(const casacore::MeasurementSet &ms,
                               const ElementResponseModel model,
                               const Options &options) {
  // Read Telescope name and convert to enum
  casacore::ScalarColumn<casacore::String> telescope_name_col(ms.observation(),
                                                              "TELESCOPE_NAME");

  TelescopeType telescope_name = Convert(telescope_name_col(0));
  switch (telescope_name) {
    case LOFAR_TELESCOPE: {
      auto telescope =
          telescope::Telescope::Ptr(new telescope::LOFAR(ms, model, options));
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