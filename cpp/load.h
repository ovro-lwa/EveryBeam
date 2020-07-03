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
#include "./telescope/lofar.h"
#include "options.h"

namespace everybeam {

/**
 * @brief Load telescope given a measurement set. Telescope is determined
 * from MeasurementSet meta-data.
 *
 * @param ms MeasurementSet
 * @param model Element response model
 * @param options Options
 * @return telescope::Telescope::Ptr
 */
telescope::Telescope::Ptr Load(const casacore::MeasurementSet &ms,
                               const ElementResponseModel model,
                               const Options &options = Options::GetDefault());
}  // namespace everybeam
#endif  // EVERYBEAM_LOAD_H_