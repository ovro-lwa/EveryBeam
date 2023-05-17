// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../cpp/elementresponse.h"

#include "config.h"

#include "tElementBeamCommon.h"

using everybeam::ElementResponseModel;

int main() {
  ElementResponseModel model(ElementResponseModel::kHamaker);
  double frequency = 132e6;  // Mhz
  std::string input_filename(TEST_MEASUREMENTSET);
  std::string output_filename("element-beams-hamaker.fits");
  run(model, frequency, input_filename, output_filename);
}