// lofarreadutils.cc: Utility functions to read the meta data relevant for
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

#include "lofarreadutils.h"
#include "beamformeridenticalantennas.h"
#include "common/mathutils.h"
#include "common/casautils.h"

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

#include <cassert>
#include <stdexcept>

#include <casacore/ms/MeasurementSets/MSAntenna.h>
#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/ms/MSSel/MSAntennaParse.h>
#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/ms/MeasurementSets/MSDataDescription.h>
#include <casacore/ms/MeasurementSets/MSDataDescColumns.h>
#include <casacore/ms/MeasurementSets/MSField.h>
#include <casacore/ms/MeasurementSets/MSFieldColumns.h>
#include <casacore/ms/MeasurementSets/MSObservation.h>
#include <casacore/ms/MeasurementSets/MSObsColumns.h>
#include <casacore/ms/MeasurementSets/MSPolarization.h>
#include <casacore/ms/MeasurementSets/MSPolColumns.h>
#include <casacore/ms/MeasurementSets/MSSpectralWindow.h>
#include <casacore/ms/MeasurementSets/MSSpWindowColumns.h>

namespace everybeam {

constexpr Antenna::CoordinateSystem::Axes antenna_orientation = {
    {
        1.0,
        0.0,
        0.0,
    },
    {
        0.0,
        1.0,
        0.0,
    },
    {0.0, 0.0, 1.0},
};

using namespace casacore;

vector3r_t TransformToFieldCoordinates(
    const vector3r_t &position, const Antenna::CoordinateSystem::Axes &axes);

// inline Antenna::CoordinateSystem ReadCoordinateSystem(
//     const casacore::Table &table, unsigned int id) {
//   casacore::ArrayQuantColumn<casacore::Double> c_position(table, "POSITION",
//                                                           "m");
//   casacore::ArrayQuantColumn<casacore::Double> c_axes(table,
//                                                       "COORDINATE_SYSTEM",
//                                                       "m");
//
//   // Read antenna field center (ITRF).
//   casacore::Vector<casacore::Quantity> aips_position = c_position(id);
//   assert(aips_position.size() == 3);
//
//   vector3r_t position = {{aips_position(0).getValue(),
//                           aips_position(1).getValue(),
//                           aips_position(2).getValue()}};
//
//   // Read antenna field coordinate axes (ITRF).
//   casacore::Matrix<casacore::Quantity> aips_axes = c_axes(id);
//   assert(aips_axes.shape().isEqual(casacore::IPosition(2, 3, 3)));
//
//   vector3r_t p = {{aips_axes(0, 0).getValue(), aips_axes(1, 0).getValue(),
//                    aips_axes(2, 0).getValue()}};
//   vector3r_t q = {{aips_axes(0, 1).getValue(), aips_axes(1, 1).getValue(),
//                    aips_axes(2, 1).getValue()}};
//   vector3r_t r = {{aips_axes(0, 2).getValue(), aips_axes(1, 2).getValue(),
//                    aips_axes(2, 2).getValue()}};
//
//   Antenna::CoordinateSystem coordinate_system = {position, {p, q, r}};
//   return coordinate_system;
// }

BeamFormer::Ptr ReadMSv3AntennaField(const Table &table, unsigned int id,
                                     ElementResponse::Ptr element_response) {
  Antenna::CoordinateSystem coordinate_system =
      common::ReadCoordinateSystem(table, id);
  BeamFormer::Ptr beam_former(
      new BeamFormerIdenticalAntennas(coordinate_system));
  //   BeamFormer::Ptr beam_former(new BeamFormer(coordinate_system));

  ROArrayQuantColumn<Double> c_offset(table, "ELEMENT_OFFSET", "m");
  ROArrayColumn<Bool> c_flag(table, "ELEMENT_FLAG");

  // Read element offsets and flags.
  Matrix<Quantity> aips_offset = c_offset(id);

  assert(aips_offset.shape().isEqual(IPosition(2, aips_offset.nrow(), 3)));

  Matrix<Bool> aips_flag = c_flag(id);
  assert(aips_flag.shape().isEqual(IPosition(2, aips_offset.nrow(), 2)));

  for (size_t i = 0; i < aips_offset.nrow(); ++i) {
    vector3r_t antenna_position = {aips_offset(i, 0).getValue(),
                                   aips_offset(i, 1).getValue(),
                                   aips_offset(i, 2).getValue()};
    antenna_position =
        TransformToFieldCoordinates(antenna_position, coordinate_system.axes);
    Antenna::Ptr antenna;
    Antenna::CoordinateSystem antenna_coordinate_system{antenna_position,
                                                        antenna_orientation};
    antenna = Element::Ptr(
        new Element(antenna_coordinate_system, element_response, id));

    antenna->enabled_[0] = !aips_flag(i, 0);
    antenna->enabled_[1] = !aips_flag(i, 1);
    beam_former->AddAntenna(antenna);
  }
  return beam_former;
}

vector3r_t ReadStationPhaseReference(const Table &table, unsigned int id);
// {
//   vector3r_t phase_reference = {0.0, 0.0, 0.0};
//   const string columnName("LOFAR_PHASE_REFERENCE");
//   if (common::HasColumn(table, columnName)) {
//     ROScalarMeasColumn<MPosition> c_reference(table, columnName);
//     MPosition mReference =
//         MPosition::Convert(c_reference(id), MPosition::ITRF)();
//     MVPosition mvReference = mReference.getValue();
//     phase_reference = {mvReference(0), mvReference(1), mvReference(2)};
//   }
//   return phase_reference;
// }

Station::Ptr ReadMSv3Station(const MeasurementSet &ms, unsigned int id,
                             const ElementResponseModel model) {
  ROMSAntennaColumns antenna(ms.antenna());
  assert(antenna.nrow() > id && !antenna.flagRow()(id));

  // Get station name.
  const string name(antenna.name()(id));

  // Get station position (ITRF).
  MPosition mPosition =
      MPosition::Convert(antenna.positionMeas()(id), MPosition::ITRF)();
  MVPosition mvPosition = mPosition.getValue();
  const vector3r_t position = {{mvPosition(0), mvPosition(1), mvPosition(2)}};

  // Create station.
  Station::Ptr station(new Station(name, position, model));

  // Read phase reference position (if available).
  //   station->SetPhaseReference(ReadStationPhaseReference(ms.antenna(), id));

  Table tab_phased_array = common::GetSubTable(ms, "PHASED_ARRAY");

  // The Station will consist of a BeamFormer that combines the fields
  // coordinate system is ITRF
  // phase reference is station position
  auto beam_former =
      ReadMSv3AntennaField(tab_phased_array, id, station->GetElementResponse());

  // TODO
  // If There is only one field, the top level beamformer is not needed
  // and the station antenna can be set the the beamformer of the field

  station->SetAntenna(beam_former);

  size_t field_id = 0;
  size_t element_id = 0;
  Antenna::CoordinateSystem coordinate_system =
      common::ReadCoordinateSystem(tab_phased_array, field_id);
  auto element_response = station->GetElementResponse();
  // TODO: rotate coordinate system for antenna
  auto element = Element::Ptr(
      new Element(coordinate_system, element_response, element_id));
  station->SetElement(element);

  return station;
}

}  // namespace everybeam
