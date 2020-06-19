// LofarMetaDataUtil.h: Utility functions to read the meta data relevant for
// simulating the beam from LOFAR observations stored in MS format.
//
// Copyright (C) 2013
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//
// $Id$

#ifndef EVERYBEAM_LOFARMETADATAUTIL_H
#define EVERYBEAM_LOFARMETADATAUTIL_H

// \file
// Utility functions to read the meta data relevant for simulating the beam from
// LOFAR observations stored in MS format.

#include "Station.h"
#include "ElementResponse.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/measures/Measures/MDirection.h>

namespace everybeam {

const ElementResponseModel defaultElementResponseModel =
    ElementResponseModel::Hamaker;

/**
 * @brief Read single station from MeasurementSet
 *
 * @param ms Measurement set
 * @param id Station id
 * @param model Element response model
 * @return Station::Ptr
 */
Station::Ptr readStation(
    const casacore::MeasurementSet &ms, unsigned int id,
    const ElementResponseModel model = defaultElementResponseModel);

/**
 * @brief Read multiple stations from measurment set into buffer out_it
 * Loops over readStation for all the antennas in MeasurementSet
 *
 * @tparam T Template type
 * @param ms Measurement set
 * @param out_it Out buffer
 * @param model Element Response buffer
 */
template <typename T>
void readStations(
    const casacore::MeasurementSet &ms, T out_it,
    const ElementResponseModel model = defaultElementResponseModel) {
  casacore::ROMSAntennaColumns antenna(ms.antenna());
  for (unsigned int i = 0; i < antenna.nrow(); ++i) {
    *out_it++ = readStation(ms, i, model);
  }
}

// Read the tile beam direction from a LOFAR MS. If it is not defined,
// this function returns the delay center.
casacore::MDirection readTileBeamDirection(const casacore::MeasurementSet &ms);

}  // namespace everybeam

#endif
