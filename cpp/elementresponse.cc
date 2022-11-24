// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// Note: ElementResponse::GetInstance is implemented in
// elementresponsefactory.cc instead of in this file. Putting it in a separate
// file removes dependencies on HamakerElementResponse::GetInstance,
// OSKARElementResponseDipole::GetInstance, LOBESElementResponse::GetInstance
// etc, which allows linking adding elementreponse.cc to libeverybeam-oskar,
// which does not link without elementreponse.cc.
// TODO (AST-1021): Find a better way of handling library dependencies.

#include "elementresponse.h"

#include <algorithm>
#include <cstdlib>

#include "config.h"
#include "options.h"

#include "common/mathutils.h"

#include "elementresponsefixeddirection.h"

namespace everybeam {

ElementResponseModel ElementResponseModelFromString(
    const std::string& element_response) {
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

std::ostream& operator<<(std::ostream& os, ElementResponseModel model) {
  switch (model) {
    case ElementResponseModel::kDefault:
      os << "Default";
      break;
    case ElementResponseModel::kHamaker:
      os << "Hamaker";
      break;
    case ElementResponseModel::kLOBES:
      os << "LOBES";
      break;
    case ElementResponseModel::kOSKARDipole:
      os << "OSKARDipole";
      break;
    case ElementResponseModel::kOSKARSphericalWave:
      os << "OSKARSphericalWave";
      break;
    case ElementResponseModel::kSkaMidAnalytical:
      os << "SKA MID Analytical Beam";
      break;
    default:
      os.setstate(std::ios_base::failbit);
  }
  return os;
}

std::shared_ptr<ElementResponse> ElementResponse::FixateDirection(
    const vector3r_t& direction) const {
  const vector2r_t thetaphi = cart2thetaphi(direction);
  return std::make_shared<ElementResponseFixedDirection>(
      shared_from_this(), thetaphi[0], thetaphi[1]);
}

std::filesystem::path ElementResponse::GetPath(
    const std::filesystem::path& relative_path) {
  return GetDataDirectory() / relative_path;
}

}  // namespace everybeam
