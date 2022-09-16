// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

// elementresponse.cc describes why ElementResponse::GetInstance is in this
// separate file.

#include "elementresponse.h"

#include <iostream>

#include "hamaker/hamakerelementresponse.h"
#include "oskar/oskarelementresponse.h"
#include "lobes/lobeselementresponse.h"
#include "options.h"

namespace everybeam {

std::shared_ptr<ElementResponse> ElementResponse::GetInstance(
    ElementResponseModel model, const std::string& name,
    const Options& options) {
  switch (model) {
    case ElementResponseModel::kHamaker:
      return HamakerElementResponse::GetInstance(name);
    case ElementResponseModel::kHamakerLba:
      return HamakerElementResponse::GetLbaInstance();
    case ElementResponseModel::kOSKARDipole:
      return OSKARElementResponseDipole::GetInstance();
    case ElementResponseModel::kOSKARSphericalWave:
      return OSKARElementResponseSphericalWave::GetInstance();
    case ElementResponseModel::kLOBES:
      try {
        return LOBESElementResponse::GetInstance(name, options);
      } catch (const std::runtime_error& e) {
        std::cout << "Creating LOBESElementResponse for station " << name
                  << " failed because: " << std::endl;
        std::cout << e.what() << std::endl;
        std::cout << "Switching to HamakerElementResponse instead" << std::endl;
        return GetInstance(ElementResponseModel::kHamaker, name, options);
      }
    default:
      std::stringstream message;
      message << "The requested element response model '" << model
              << "' is not implemented.";
      throw std::runtime_error(message.str());
  }
}

}  // namespace everybeam
