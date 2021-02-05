// load_telescope.h: Main interface function for loading a telescope
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_LOAD_H_
#define EVERYBEAM_LOAD_H_

#include "telescope/telescope.h"
#include "options.h"
#include "elementresponse.h"

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

/**
 * @brief Convert a string to an ElementResponseModel enum
 *
 */
everybeam::ElementResponseModel GetElementResponseEnum(
    const std::string &element_response);

}  // namespace everybeam

#endif  // EVERYBEAM_LOAD_H_
