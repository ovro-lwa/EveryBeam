// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "elementresponse.h"
#include "options.h"

#include "hamaker/hamakerelementresponse.h"
#include "oskar/oskarelementresponse.h"
#include "lobes/lobeselementresponse.h"

namespace everybeam {

ElementResponseModel ElementResponseModelFromString(
    const std::string &element_response) {
  // LOFAR & SKA(/OSKAR) related
  std::string element_response_upper = element_response;
  std::transform(element_response_upper.begin(), element_response_upper.end(),
                 element_response_upper.begin(), ::toupper);

  everybeam::ElementResponseModel element_response_enum;
  if (element_response_upper == "" || element_response_upper == "DEFAULT") {
    element_response_enum = everybeam::ElementResponseModel::kDefault;
  } else if (element_response_upper == "HAMAKER") {
    element_response_enum = everybeam::ElementResponseModel::kHamaker;
  } else if (element_response_upper == "LOBES") {
    element_response_enum = everybeam::ElementResponseModel::kLOBES;
  } else if (element_response_upper == "OSKARDIPOLE") {
    element_response_enum = everybeam::ElementResponseModel::kOSKARDipole;
  } else if (element_response_upper == "OSKARSPHERICALWAVE") {
    element_response_enum =
        everybeam::ElementResponseModel::kOSKARSphericalWave;
  } else {
    std::stringstream message;
    message << "The specified element response model " << element_response
            << " is not implemented.";
    throw std::runtime_error(message.str());
  }
  return element_response_enum;
}

std::ostream &operator<<(std::ostream &os, ElementResponseModel model) {
  switch (model) {
    case kDefault:
      os << "Default";
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

std::shared_ptr<ElementResponse> ElementResponse::GetInstance(
    ElementResponseModel model, const std::string &name,
    const Options &options) {
  switch (model) {
    case kHamaker:
      return HamakerElementResponse::GetInstance(name);
    case kOSKARDipole:
      return OSKARElementResponseDipole::GetInstance();
    case kOSKARSphericalWave:
      return OSKARElementResponseSphericalWave::GetInstance();
    case kLOBES:
      try {
        return LOBESElementResponse::GetInstance(name, options);
      } catch (const std::runtime_error &e) {
        std::cout << "Creating LOBESElementResponse for station " << name
                  << " failed because: " << std::endl;
        std::cout << e.what() << std::endl;
        std::cout << "Switching to HamakerElementResponse instead" << std::endl;
        return GetInstance(kHamaker, name, options);
      }
    default:
      std::stringstream message;
      message << "The requested element response model '" << model
              << "' is not implemented.";
      throw std::runtime_error(message.str());
  }
}

}  // namespace everybeam
