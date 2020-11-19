// lofarreadutils.h: Utility functions to read the meta data relevant for
// simulating the beam from LOFAR observations stored in MS format.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_LOFARREADUTILS_H_
#define EVERYBEAM_LOFARREADUTILS_H_

// \file
// Utility functions to read the meta data relevant for simulating the beam from
// LOFAR observations stored in MS format.

#include "station.h"
#include "elementresponse.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/measures/Measures/MDirection.h>

namespace everybeam {
const ElementResponseModel defaultElementResponseModel =
    ElementResponseModel::kHamaker;

/**
 * @brief Read single station from MeasurementSet
 *
 * @param ms Measurement set
 * @param id Station id
 * @param model Element response model
 * @return Station::Ptr
 */
Station::Ptr ReadLofarStation(
    const casacore::MeasurementSet &ms, unsigned int id,
    const ElementResponseModel model = defaultElementResponseModel);

/**
 * @brief Read multiple stations from measurment set into buffer out_it
 * Loops over ReadLofarStation for all the antennas in MeasurementSet
 *
 * @tparam T Template type
 * @param ms Measurement set
 * @param out_it Out buffer
 * @param model Element Response buffer
 */
template <typename T>
void ReadStations(
    const casacore::MeasurementSet &ms, T out_it,
    const ElementResponseModel model = defaultElementResponseModel) {
  casacore::ROMSAntennaColumns antenna(ms.antenna());
  for (unsigned int i = 0; i < antenna.nrow(); ++i) {
    *out_it++ = ReadLofarStation(ms, i, model);
  }
}

// Read the tile beam direction from a LOFAR MS. If it is not defined,
// this function returns the delay center.
casacore::MDirection ReadTileBeamDirection(const casacore::MeasurementSet &ms);
}  // namespace everybeam
#endif  // EVERYBEAM_LOFARREADUTILS_H_
