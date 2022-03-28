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
/** AARTFAAC antenna column name in MS */
const std::string kAartfaacAntennaTypeName = "AARTFAAC_ANTENNA_TYPE";

/**
 * @brief Read single station from MeasurementSet by index
 *
 * @param ms Measurement set
 * @param id Station id
 * @param model Element response model
 * @param options [optional] can contain for example the coefficient path, used
 to specify
 * locations of LOBES or other coefficient file(s)
 * @return Shared pointer to \param Station object

 */
std::shared_ptr<Station> ReadSingleStation(const casacore::MeasurementSet& ms,
                                           unsigned int id,
                                           const Options& options = Options());

/**
 * @brief Read multiple stations from measurment set into buffer out_it
 * Loops over ReadSingleStation for all the antennas in MeasurementSet
 *
 * @tparam T Template type
 * @param ms Measurement set
 * @param out_it Out buffer, storing shared pointers to Station objects
 * @param model ElementResponseModel to use
 */
template <typename T>
inline void ReadAllStations(const casacore::MeasurementSet& ms, T out_it,
                            const Options& options = Options()) {
  casacore::ROMSAntennaColumns antenna(ms.antenna());
  for (unsigned int i = 0; i < antenna.nrow(); ++i) {
    *out_it++ = ReadSingleStation(ms, i, options);
  }
}

// Read the tile beam direction from a LOFAR MS. If it is not defined,
// this function returns the delay center.
casacore::MDirection ReadTileBeamDirection(const casacore::MeasurementSet& ms);
}  // namespace everybeam
#endif  // EVERYBEAM_MSREADUTILS_H_