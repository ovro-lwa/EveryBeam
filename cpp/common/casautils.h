// casautils.h: CasaCore utilities.
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

#ifndef EVERYBEAM_COMMON_CASAUTIL_H_
#define EVERYBEAM_COMMON_CASAUTIL_H_

#include "types.h"
#include "./../antenna.h"

#include <cassert>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/tables/Tables/TableRecord.h>

namespace everybeam {
namespace common {

/**
 * @brief Read coordinate system from MeasurementSet
 *
 * @param table Measurement set (casacore::Table)
 * @param id Id of the antenna field in the station (int)
 * @return Antenna::CoordinateSystem
 */
inline Antenna::CoordinateSystem ReadCoordinateSystem(
    const casacore::Table &table, unsigned int id) {
  casacore::ArrayQuantColumn<casacore::Double> c_position(table, "POSITION",
                                                          "m");
  casacore::ArrayQuantColumn<casacore::Double> c_axes(table, "COORDINATE_AXES",
                                                      "m");

  // Read antenna field center (ITRF).
  casacore::Vector<casacore::Quantity> aips_position = c_position(id);
  assert(aips_position.size() == 3);

  vector3r_t position = {{aips_position(0).getValue(),
                          aips_position(1).getValue(),
                          aips_position(2).getValue()}};

  // Read antenna field coordinate axes (ITRF).
  casacore::Matrix<casacore::Quantity> aips_axes = c_axes(id);
  assert(aips_axes.shape().isEqual(casacore::IPosition(2, 3, 3)));

  vector3r_t p = {{aips_axes(0, 0).getValue(), aips_axes(1, 0).getValue(),
                   aips_axes(2, 0).getValue()}};
  vector3r_t q = {{aips_axes(0, 1).getValue(), aips_axes(1, 1).getValue(),
                   aips_axes(2, 1).getValue()}};
  vector3r_t r = {{aips_axes(0, 2).getValue(), aips_axes(1, 2).getValue(),
                   aips_axes(2, 2).getValue()}};

  Antenna::CoordinateSystem coordinate_system = {position, {p, q, r}};
  return coordinate_system;
}

/**
 * @brief Check if the specified column exists as a column of the
 * specified table.
 *
 * @param table Measurement set (casacore::Table)
 * @param column Column name (str)
 * @return true If column present
 * @return false If column not present
 */
inline bool HasColumn(const casacore::Table &table, const string &column) {
  return table.tableDesc().isColumn(column);
}

/**
 * @brief Provide access to a sub-table by name.
 *
 * @param table Measurment set (casacore::Table)
 * @param name Name of sub table (str)
 * @return Table (casacore::Table)
 */
inline casacore::Table GetSubTable(const casacore::Table &table,
                                   const string &name) {
  return table.keywordSet().asTable(name);
}
}  // namespace common
}  // namespace everybeam
#endif  // EVERYBEAM_COMMON_CASAUTIL_H_