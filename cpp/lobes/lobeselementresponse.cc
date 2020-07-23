#include "lobeselementresponse.h"

#include <map>

namespace everybeam {

LOBESElementResponse::LOBESElementResponse(std::string name) {}

void LOBESElementResponse::Response(
    double freq, double theta, double phi,
    std::complex<double> (&response)[2][2]) const {}

std::shared_ptr<LOBESElementResponse> LOBESElementResponse::GetInstance(
    std::string name) {
  static std::map<std::string, std::shared_ptr<LOBESElementResponse>>
      name_response_map;

  auto entry = name_response_map.find(name);
  if (entry == name_response_map.end()) {
    entry = name_response_map.insert(
        entry, {name, std::make_shared<LOBESElementResponse>(name)});
  }
  return entry->second;
};

}  // namespace everybeam
