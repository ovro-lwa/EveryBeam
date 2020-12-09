// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_MSREADUTILS_H_
#define EVERYBEAM_MSREADUTILS_H_

// \file
// Utility functions to read the meta data relevant for simulating the beam from
// LOFAR / OSKAR observations stored in MS format.

#include "station.h"
#include "elementresponse.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/measures/Measures/MDirection.h>

namespace everybeam {
/**
 * @brief Read single station from MeasurementSet by index
 *
 * @param ms Measurement set
 * @param id Station id
 * @param model Element response model
 * @return shared
 */
std::shared_ptr<Station> ReadSingleStation(
    const casacore::MeasurementSet &ms, unsigned int id,
    ElementResponseModel model = ElementResponseModel::kDefault);

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
inline void ReadAllStations(
    const casacore::MeasurementSet &ms, T out_it,
    const ElementResponseModel model = ElementResponseModel::kDefault) {
  casacore::ROMSAntennaColumns antenna(ms.antenna());
  for (unsigned int i = 0; i < antenna.nrow(); ++i) {
    *out_it++ = ReadSingleStation(ms, i, model);
  }
}

// Read the tile beam direction from a LOFAR MS. If it is not defined,
// this function returns the delay center.
casacore::MDirection ReadTileBeamDirection(const casacore::MeasurementSet &ms);
}  // namespace everybeam
#endif  // EVERYBEAM_MSREADUTILS_H_