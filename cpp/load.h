// load_telescope.h: Main interface function for loading a telescope
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

#ifndef EVERYBEAM_LOAD_H_
#define EVERYBEAM_LOAD_H_

#include "./telescope/telescope.h"
#include "options.h"

namespace everybeam {
/**
 * @brief Available TelescopeType enums
 *
 */
enum TelescopeType {
  kUnknownTelescope,
  kLofarTelescope,
  kAARTFAAC,
  kVLATelescope,
  kATCATelescope,
  kMWATelescope,
  kOSKARTelescope
};

/**
 * @brief Derive the TelescopeType from a given MS
 *
 * @param ms
 * @return TelescopeType
 */
TelescopeType GetTelescopeType(const casacore::MeasurementSet &ms);

/**
 * @brief Load telescope given a measurement set. Telescope is determined
 * from MeasurementSet meta-data.
 *
 * @param ms MeasurementSet
 * @param options Options
 * @return telescope::Telescope::Ptr
 */
std::unique_ptr<telescope::Telescope> Load(const casacore::MeasurementSet &ms,
                                           const Options &options);

/**
 * @brief Load telescope given a path to a measurment set. Telescope is
 * determined from MeasurementSet meta-data.
 *
 * @param ms MeasurementSet
 * @param options Options
 * @return telescope::Telescope::Ptr
 */
std::unique_ptr<telescope::Telescope> Load(const std::string &ms_name,
                                           const Options &options);
}  // namespace everybeam

#endif  // EVERYBEAM_LOAD_H_
