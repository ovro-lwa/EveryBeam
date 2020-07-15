#include "element_response.h"

namespace everybeam {
std::ostream& operator<<(std::ostream& os, ElementResponseModel model) {
  switch (model) {
    case kUnknown:
      os << "Unknown";
      break;
    case kHamaker:
      os << "Hamaker";
      break;
    case kLOBES:
      os << "LOBES";
      break;
    case kOSKARDipole:
      os << "OSKARDipole";
      break;
    case kOSKARSphericalWave:
      os << "OSKARSphericalWave";
      break;
    default:
      os.setstate(std::ios_base::failbit);
  }
  return os;
}
}  // namespace everybeam
