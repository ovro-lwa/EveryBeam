// Options.h: Class for specifying the telescope response options.
//
// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the EveryBeam software suite.
// The EveryBeam software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The EveryBeam software suite is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the EveryBeam software suite. If not, see
// <http://www.gnu.org/licenses/>.
//
// $Id$

#ifndef EVERYBEAM_OPTIONS_H_
#define EVERYBEAM_OPTIONS_H_

#include <string>
#include "elementresponse.h"

namespace everybeam {

/**
 * @brief Class/Struct specifying everybeam Options. Needs further
 * implementation!
 *
 */
struct Options {
  // Path to coefficients file
  std::string coeff_path = ".";

  // LOFAR specific
  bool use_differential_beam = false;
  bool use_channel_frequency = true;
  std::string data_column_name = "DATA";
  ElementResponseModel element_response_model = ElementResponseModel::kHamaker;

  // MWA specific (Lofar probably will follow)
  bool frequency_interpolation = false;
};
}  // namespace everybeam
#endif  // EVERYBEAM_OPTIONS_H_