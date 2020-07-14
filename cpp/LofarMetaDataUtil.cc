// LofarMetaDataUtil.cc: Utility functions to read the meta data relevant for
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

#include "LofarMetaDataUtil.h"
#include "common/math_utils.h"
#include "common/casa_utils.h"

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

constexpr Antenna::CoordinateSystem::Axes lofar_antenna_orientation = {
    {
        -std::sqrt(.5),
        -std::sqrt(.5),
        0.0,
    },
    {
        std::sqrt(.5),
        -std::sqrt(.5),
        0.0,
    },
    {0.0, 0.0, 1.0},
};

// constexpr Antenna::CoordinateSystem::Axes lofar_antenna_orientation = {
//     {1.0, 0.0, 0.0,},
//     {0.0, 1.0, 0.0,},
//     {0.0, 0.0, 1.0},
// };

using namespace casacore;

typedef std::array<vector3r_t, 16> TileConfig;

// TODO: utility is not used at all
bool hasSubTable(const Table &table, const string &name) {
  return table.keywordSet().isDefined(name);
}

TileConfig readTileConfig(const Table &table, unsigned int row) {
  ROArrayQuantColumn<Double> c_tile_offset(table, "TILE_ELEMENT_OFFSET", "m");

  // Read tile configuration for HBA antenna fields, assert validity of aips
  // offset.
  Matrix<Quantity> aips_offset = c_tile_offset(row);

  TileConfig config;
  assert(aips_offset.ncolumn() == config.size());

  for (unsigned int i = 0; i < config.size(); ++i) {
    config[i][0] = aips_offset(0, i).getValue();
    config[i][1] = aips_offset(1, i).getValue();
    config[i][2] = aips_offset(2, i).getValue();
  }
  return config;
}

void transformToFieldCoordinates(TileConfig &config,
                                 const Antenna::CoordinateSystem::Axes &axes) {
  for (unsigned int i = 0; i < config.size(); ++i) {
    const vector3r_t position = config[i];
    config[i][0] = dot(position, axes.p);
    config[i][1] = dot(position, axes.q);
    config[i][2] = dot(position, axes.r);
  }
}

vector3r_t transformToFieldCoordinates(
    const vector3r_t &position, const Antenna::CoordinateSystem::Axes &axes) {
  const vector3r_t result{dot(position, axes.p), dot(position, axes.q),
                          dot(position, axes.r)};
  return result;
}

// AntennaField::CoordinateSystem readCoordinateSystemAartfaac(
//     const Table &table, unsigned int id)
// {
//     ROArrayQuantColumn<Double> c_position(table, "POSITION", "m");
//
//     // Read antenna field center (ITRF).
//     Vector<Quantity> aips_position = c_position(id);
//     assert(aips_position.size() == 3);
//
//     vector3r_t position = {{aips_position(0).getValue(),
//         aips_position(1).getValue(), aips_position(2).getValue()}};
//
//     TableRecord keywordset = table.keywordSet();
//     Matrix<double> aips_axes;
//     keywordset.get("AARTFAAC_COORDINATE_AXES", aips_axes);
//     assert(aips_axes.shape().isEqual(IPosition(2, 3, 3)));
//
//     vector3r_t p = {{aips_axes(0, 0), aips_axes(1, 0), aips_axes(2, 0)}};
//     vector3r_t q = {{aips_axes(0, 1), aips_axes(1, 1), aips_axes(2, 1)}};
//     vector3r_t r = {{aips_axes(0, 2), aips_axes(1, 2), aips_axes(2, 2)}};
//
//     AntennaField::CoordinateSystem system = {position, {p, q, r}};
//
//     return system;
// }

Antenna::CoordinateSystem readCoordinateSystem(const Table &table,
                                               unsigned int id) {
  ROArrayQuantColumn<Double> c_position(table, "POSITION", "m");
  ROArrayQuantColumn<Double> c_axes(table, "COORDINATE_AXES", "m");

  // Read antenna field center (ITRF).
  Vector<Quantity> aips_position = c_position(id);
  assert(aips_position.size() == 3);

  vector3r_t position = {{aips_position(0).getValue(),
                          aips_position(1).getValue(),
                          aips_position(2).getValue()}};

  // Read antenna field coordinate axes (ITRF).
  Matrix<Quantity> aips_axes = c_axes(id);
  assert(aips_axes.shape().isEqual(IPosition(2, 3, 3)));

  vector3r_t p = {{aips_axes(0, 0).getValue(), aips_axes(1, 0).getValue(),
                   aips_axes(2, 0).getValue()}};
  vector3r_t q = {{aips_axes(0, 1).getValue(), aips_axes(1, 1).getValue(),
                   aips_axes(2, 1).getValue()}};
  vector3r_t r = {{aips_axes(0, 2).getValue(), aips_axes(1, 2).getValue(),
                   aips_axes(2, 2).getValue()}};

  Antenna::CoordinateSystem coordinate_system = {position, {p, q, r}};
  return coordinate_system;
}

BeamFormer::Ptr make_tile(unsigned int id, const vector3r_t &position,
                          const TileConfig &tile_config,
                          ElementResponse::Ptr element_response) {
  BeamFormer::Ptr tile = BeamFormer::Ptr(new BeamFormer(position));

  for (unsigned int id = 0; id < tile_config.size(); id++) {
    vector3r_t antenna_position = tile_config[id];

    Antenna::CoordinateSystem antenna_coordinate_system;
    antenna_coordinate_system.origin = antenna_position;
    antenna_coordinate_system.axes = lofar_antenna_orientation;

    Antenna::Ptr antenna = Element::Ptr(
        new Element(antenna_coordinate_system, element_response, id));
    tile->AddAntenna(antenna);
  }

  return tile;
}

BeamFormer::Ptr readAntennaField(const Table &table, unsigned int id,
                                 ElementResponse::Ptr element_response) {
  Antenna::CoordinateSystem coordinate_system =
      common::readCoordinateSystem(table, id);
  //     std::cout << "coordinate_system: " << std::endl;
  //     std::cout << "  axes.p: " << coordinate_system.axes.p[0] << ", " <<
  //     coordinate_system.axes.p[1] << ", " << coordinate_system.axes.p[2] <<
  //     std::endl; std::cout << "  axes.q: " << coordinate_system.axes.q[0] <<
  //     ", " << coordinate_system.axes.q[1] << ", " <<
  //     coordinate_system.axes.q[2] << std::endl; std::cout << "  axes.r: " <<
  //     coordinate_system.axes.r[0] << ", " << coordinate_system.axes.r[1] <<
  //     ", " << coordinate_system.axes.r[2] << std::endl;
  BeamFormer::Ptr beam_former(new BeamFormer(coordinate_system));

  ROScalarColumn<String> c_name(table, "NAME");
  ROArrayQuantColumn<Double> c_offset(table, "ELEMENT_OFFSET", "m");
  ROArrayColumn<Bool> c_flag(table, "ELEMENT_FLAG");

  const string &name = c_name(id);

  // Read element offsets and flags.
  Matrix<Quantity> aips_offset = c_offset(id);
  assert(aips_offset.shape().isEqual(IPosition(2, 3, aips_offset.ncolumn())));

  Matrix<Bool> aips_flag = c_flag(id);
  assert(aips_flag.shape().isEqual(IPosition(2, 2, aips_offset.ncolumn())));

  TileConfig tile_config;
  if (name != "LBA") tile_config = readTileConfig(table, id);
  transformToFieldCoordinates(tile_config, coordinate_system.axes);

  for (size_t i = 0; i < aips_offset.ncolumn(); ++i) {
    vector3r_t antenna_position = {aips_offset(0, i).getValue(),
                                   aips_offset(1, i).getValue(),
                                   aips_offset(2, i).getValue()};
    antenna_position =
        transformToFieldCoordinates(antenna_position, coordinate_system.axes);
    Antenna::Ptr antenna;
    Antenna::CoordinateSystem antenna_coordinate_system{
        antenna_position, lofar_antenna_orientation};
    if (name == "LBA") {
      antenna = Element::Ptr(
          new Element(antenna_coordinate_system, element_response, id));
    } else {
      // name is HBA, HBA0, HBA1
      antenna = make_tile(id, antenna_position, tile_config, element_response);
    }

    antenna->m_enabled[0] = !aips_flag(0, i);
    antenna->m_enabled[1] = !aips_flag(1, i);
    beam_former->AddAntenna(antenna);
  }
  return beam_former;
}

BeamFormer::Ptr readAntennaFieldAartfaac(const Table &table,
                                         const string &ant_type,
                                         unsigned int id) {
  BeamFormer::Ptr field;
  //     AntennaField::CoordinateSystem system =
  //     readCoordinateSystemAartfaac(table, id);
  //
  //     if (ant_type == "LBA")
  //     {
  //         DualDipoleAntenna::Ptr model(new DualDipoleAntenna());
  //         field = AntennaField::Ptr(new AntennaFieldLBA(ant_type, system,
  //         model));
  //     }
  //     else // HBA
  //     {
  //          // TODO: implement this
  //          throw std::runtime_error("HBA for Aartfaac is not implemented
  //          yet.");
  //     }
  //
  //     // Add only one antenna to the field (no offset, always enabled)
  //     AntennaField::Antenna antenna;
  //     antenna.position[0] = 0.;
  //     antenna.position[1] = 0.;
  //     antenna.position[2] = 0.;
  //     antenna.enabled[0] = true;
  //     antenna.enabled[1] = true;
  //
  //     field->addAntenna(antenna);

  return field;
}

vector3r_t readStationPhaseReference(const Table &table, unsigned int id) {
  vector3r_t phase_reference = {0.0, 0.0, 0.0};
  const string columnName("LOFAR_PHASE_REFERENCE");
  if (common::hasColumn(table, columnName)) {
    ROScalarMeasColumn<MPosition> c_reference(table, columnName);
    MPosition mReference =
        MPosition::Convert(c_reference(id), MPosition::ITRF)();
    MVPosition mvReference = mReference.getValue();
    phase_reference = {mvReference(0), mvReference(1), mvReference(2)};
  }
  return phase_reference;
}

Station::Ptr readStation(const MeasurementSet &ms, unsigned int id,
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
  station->setPhaseReference(readStationPhaseReference(ms.antenna(), id));

  // Read antenna field information.
  ROScalarColumn<String> telescope_name_col(
      common::getSubTable(ms, "OBSERVATION"), "TELESCOPE_NAME");
  string telescope_name = telescope_name_col(0);

  if (telescope_name == "LOFAR") {
    Table tab_field = common::getSubTable(ms, "LOFAR_ANTENNA_FIELD");
    tab_field = tab_field(tab_field.col("ANTENNA_ID") == static_cast<Int>(id));

    // The Station will consist of a BeamFormer that combines the fields
    // coordinate system is ITRF
    // phase reference is station position
    auto beam_former = BeamFormer::Ptr(new BeamFormer(
        Antenna::IdentityCoordinateSystem, station->phaseReference()));

    for (size_t i = 0; i < tab_field.nrow(); ++i) {
      beam_former->AddAntenna(
          readAntennaField(tab_field, i, station->get_element_response()));
    }

    // TODO
    // If There is only one field, the top level beamformer is not needed
    // and the station antenna can be set the the beamformer of the field

    station->SetAntenna(beam_former);

    size_t field_id = 0;
    size_t element_id = 0;
    Antenna::CoordinateSystem coordinate_system =
        common::readCoordinateSystem(tab_field, field_id);
    auto model = station->get_element_response();
    // TODO: rotate coordinate system for antenna
    auto element =
        Element::Ptr(new Element(coordinate_system, model, element_id));
    station->SetElement(element);
  } else if (telescope_name == "AARTFAAC") {
    ROScalarColumn<String> ant_type_col(common::getSubTable(ms, "OBSERVATION"),
                                        "AARTFAAC_ANTENNA_TYPE");
    string ant_type = ant_type_col(0);

    Table tab_field = common::getSubTable(ms, "ANTENNA");
    station->SetAntenna(readAntennaFieldAartfaac(tab_field, ant_type, id));
  }

  return station;
}

MDirection readTileBeamDirection(const casacore::MeasurementSet &ms) {
  MDirection tileBeamDir;

  Table fieldTable = common::getSubTable(ms, "FIELD");

  if (fieldTable.nrow() != 1) {
    throw std::runtime_error(
        "MS has multiple fields, this does not work with the LOFAR beam "
        "library.");
  }

  if (common::hasColumn(fieldTable, "LOFAR_TILE_BEAM_DIR")) {
    ROArrayMeasColumn<MDirection> tileBeamCol(fieldTable,
                                              "LOFAR_TILE_BEAM_DIR");
    tileBeamDir = *(tileBeamCol(0).data());
  } else {
    ROArrayMeasColumn<MDirection> tileBeamCol(fieldTable, "DELAY_DIR");
    tileBeamDir = *(tileBeamCol(0).data());
  }

  return tileBeamDir;
}

}  // namespace everybeam
