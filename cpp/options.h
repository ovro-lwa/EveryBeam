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

namespace everybeam {

/**
 * @brief Class/Struct specifying everybeam Options. Needs further
 * implementation!
 *
 */
class Options {
 public:
  Options()
      : use_differential_beam(false),
        use_channel_frequency(true),
        data_column_name("DATA"){};

  //! Default - empty - options class
  static Options GetDefault() { return Options(); };

  // TODO? Specify path to element response coefficients file
  // std::string coeff_path;
  bool use_differential_beam;
  bool use_channel_frequency;
  std::string data_column_name;
};
}  // namespace everybeam
#endif  // EVERYBEAM_OPTIONS_H_