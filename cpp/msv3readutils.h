// msv3readutils.h: Utility functions to read the meta data relevant for
// simulating the beam from OSKAR simulations stored in MS format.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_MSV3READUTILS_H_
#define EVERYBEAM_MSV3READUTILS_H_

// \file
// Utility functions to read the meta data relevant for simulating the beam from
// OSKAR simulations stored in MS format.

#include "station.h"
#include "elementresponse.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/measures/Measures/MDirection.h>

namespace everybeam {
const ElementResponseModel defaultElementResponseModel =
    ElementResponseModel::kDefault;

/**
 * @brief Read single station from MeasurementSet
 *
 * @param ms Measurement set
 * @param id Station id
 * @param model Element response model
 * @return Station::Ptr
 */
Station::Ptr ReadMSv3Station(
    const casacore::MeasurementSet &ms, unsigned int id,
    const ElementResponseModel model = defaultElementResponseModel);

/**
 * @brief Read multiple stations from measurment set into buffer out_it
 * Loops over ReadMSv3Station for all the antennas in MeasurementSet
 *
 * @tparam T Template type
 * @param ms Measurement set
 * @param out_it Out buffer
 * @param model Element Response buffer
 */
template <typename T>
void ReadMSv3Stations(
    const casacore::MeasurementSet &ms, T out_it,
    const ElementResponseModel model = defaultElementResponseModel) {
  casacore::ROMSAntennaColumns antenna(ms.antenna());
  for (unsigned int i = 0; i < antenna.nrow(); ++i) {
    *out_it++ = ReadMSv3Station(ms, i, model);
  }
}

}  // namespace everybeam
#endif  // EVERYBEAM_MSV3READUTILS_H_
