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
#include "beamformerlofarhba.h"
#include "common/mathutils.h"
#include "common/casautils.h"

#include <memory>

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
bool HasSubTable(const Table &table, const string &name) {
  return table.keywordSet().isDefined(name);
}

TileConfig ReadTileConfig(const Table &table, unsigned int row) {
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

void TransformToFieldCoordinates(TileConfig &config,
                                 const Antenna::CoordinateSystem::Axes &axes) {
  for (unsigned int i = 0; i < config.size(); ++i) {
    const vector3r_t position = config[i];
    config[i][0] = dot(position, axes.p);
    config[i][1] = dot(position, axes.q);
    config[i][2] = dot(position, axes.r);
  }
}

vector3r_t TransformToFieldCoordinates(
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

std::shared_ptr<BeamFormer> MakeTile(unsigned int id,
                                     const vector3r_t &position,
                                     const TileConfig &tile_config,
                                     ElementResponse::Ptr element_response) {
  std::shared_ptr<BeamFormer> tile =
      std::make_shared<BeamFormer>(BeamFormerIdenticalAntennas(position));

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

// Make a dedicated HBA "Hamaker" tile, saving only one element, and 16
// element positions
void MakeTile(std::shared_ptr<BeamFormerLofarHBA> beamformer,
              const vector3r_t &position, const TileConfig &tile_config,
              ElementResponse::Ptr element_response) {
  for (unsigned int id = 0; id < tile_config.size(); id++) {
    vector3r_t antenna_position = tile_config[id];

    Antenna::CoordinateSystem antenna_coordinate_system;
    antenna_coordinate_system.origin = antenna_position;
    antenna_coordinate_system.axes = lofar_antenna_orientation;

    std::shared_ptr<ElementHamaker> antenna = std::make_shared<ElementHamaker>(
        ElementHamaker(antenna_coordinate_system, element_response, id));

    // Only element 1 needs to be stored as an element
    if (id == 0) {
      beamformer->SetElement(antenna);
    }
    // All positions need to be stored, however
    beamformer->AddElementPosition(antenna_position);
  }
}

Antenna::Ptr ReadAntennaField(const Table &table, unsigned int id,
                              ElementResponse::Ptr element_response,
                              ElementResponseModel element_response_model) {
  Antenna::CoordinateSystem coordinate_system =
      common::ReadCoordinateSystem(table, id);

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
  if (name != "LBA") tile_config = ReadTileConfig(table, id);

  std::shared_ptr<Antenna> beam_former;
  // Cast to the beam_former corresponding to the element response
  // model and LBA/HBA configuration
  if (element_response_model == ElementResponseModel::kHamaker) {
    if (name != "LBA") {
      // Then HBA, HBA0 or HBA1
      beam_former = std::make_shared<BeamFormerLofarHBA>(
          BeamFormerLofarHBA(coordinate_system));
    } else {
      // TODO: will become a dedicated BeamFormerLofarLBA in the near future
      beam_former = std::make_shared<BeamFormerIdenticalAntennas>(
          BeamFormerIdenticalAntennas(coordinate_system));
    }
  } else {
    // Tiles / element should be kept unique, so work with generic BeamFormer
    beam_former = std::make_shared<BeamFormer>(BeamFormer(coordinate_system));
  }

  TransformToFieldCoordinates(tile_config, coordinate_system.axes);

  for (size_t i = 0; i < aips_offset.ncolumn(); ++i) {
    vector3r_t antenna_position = {aips_offset(0, i).getValue(),
                                   aips_offset(1, i).getValue(),
                                   aips_offset(2, i).getValue()};
    antenna_position =
        TransformToFieldCoordinates(antenna_position, coordinate_system.axes);
    Antenna::Ptr antenna;
    Antenna::CoordinateSystem antenna_coordinate_system{
        antenna_position, lofar_antenna_orientation};

    if (name == "LBA") {
      antenna = std::make_shared<Element>(
          Element(antenna_coordinate_system, element_response, id));
      antenna->enabled_[0] = !aips_flag(0, i);
      antenna->enabled_[1] = !aips_flag(1, i);

      if (element_response_model == kHamaker) {
        // NOTE: no cast needed as yet, but it already hints towards
        // future implementation
        std::shared_ptr<BeamFormerIdenticalAntennas> beam_former_lba =
            std::static_pointer_cast<BeamFormerIdenticalAntennas>(beam_former);
        beam_former_lba->AddAntenna(antenna);
      } else {
        std::shared_ptr<BeamFormer> beam_former_lba =
            std::static_pointer_cast<BeamFormer>(beam_former);
        beam_former_lba->AddAntenna(antenna);
      }
    } else {
      // name is HBA, HBA0 or HBA1
      if (element_response_model == kHamaker) {
        std::shared_ptr<BeamFormerLofarHBA> beam_former_hba =
            std::static_pointer_cast<BeamFormerLofarHBA>(beam_former);

        // Tile positions are uniques
        beam_former_hba->AddTilePosition(antenna_position);

        // Store only one tile
        if (i == 0) {
          MakeTile(beam_former_hba, antenna_position, tile_config,
                   element_response);
        }
        // Tile enabled in x/y?
        beam_former_hba->AddTileEnabled(
            std::array<bool, 2>{!aips_flag(0, i), !aips_flag(0, i)});
      } else {
        std::shared_ptr<BeamFormerIdenticalAntennas> beam_former_hba =
            std::static_pointer_cast<BeamFormerIdenticalAntennas>(beam_former);
        antenna = MakeTile(id, antenna_position, tile_config, element_response);
        antenna->enabled_[0] = !aips_flag(0, i);
        antenna->enabled_[1] = !aips_flag(1, i);
        beam_former_hba->AddAntenna(antenna);
      }
    }
  }
  return beam_former;
}

BeamFormer::Ptr ReadAntennaFieldAartfaac(const Table &table,
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

vector3r_t ReadStationPhaseReference(const Table &table, unsigned int id) {
  vector3r_t phase_reference = {0.0, 0.0, 0.0};
  const string columnName("LOFAR_PHASE_REFERENCE");
  if (common::HasColumn(table, columnName)) {
    ROScalarMeasColumn<MPosition> c_reference(table, columnName);
    MPosition mReference =
        MPosition::Convert(c_reference(id), MPosition::ITRF)();
    MVPosition mvReference = mReference.getValue();
    phase_reference = {mvReference(0), mvReference(1), mvReference(2)};
  }
  return phase_reference;
}

Station::Ptr ReadLofarStation(const MeasurementSet &ms, unsigned int id,
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
  Station::Ptr station =
      std::make_shared<Station>(Station(name, position, model));

  // Read phase reference position (if available).
  station->SetPhaseReference(ReadStationPhaseReference(ms.antenna(), id));

  // Read antenna field information.
  ROScalarColumn<String> telescope_name_col(
      common::GetSubTable(ms, "OBSERVATION"), "TELESCOPE_NAME");
  string telescope_name = telescope_name_col(0);

  if (telescope_name == "LOFAR") {
    Table tab_field = common::GetSubTable(ms, "LOFAR_ANTENNA_FIELD");
    tab_field = tab_field(tab_field.col("ANTENNA_ID") == static_cast<Int>(id));

    // The Station will consist of a BeamFormer that combines the fields
    // coordinate system is ITRF
    // phase reference is station position
    auto beam_former = std::make_shared<BeamFormer>(BeamFormer(
        Antenna::IdentityCoordinateSystem, station->GetPhaseReference()));

    for (size_t i = 0; i < tab_field.nrow(); ++i) {
      beam_former->AddAntenna(
          ReadAntennaField(tab_field, i, station->GetElementResponse(),
                           station->GetElementResponseModel()));
    }

    // TODO
    // If There is only one field, the top level beamformer is not needed
    // and the station antenna can be set the the beamformer of the field
    station->SetAntenna(beam_former);
  } else if (telescope_name == "AARTFAAC") {
    ROScalarColumn<String> ant_type_col(common::GetSubTable(ms, "OBSERVATION"),
                                        "AARTFAAC_ANTENNA_TYPE");
    string ant_type = ant_type_col(0);

    Table tab_field = common::GetSubTable(ms, "ANTENNA");
    station->SetAntenna(ReadAntennaFieldAartfaac(tab_field, ant_type, id));
  }

  return station;
}

MDirection ReadTileBeamDirection(const casacore::MeasurementSet &ms) {
  MDirection tileBeamDir;

  Table fieldTable = common::GetSubTable(ms, "FIELD");

  if (fieldTable.nrow() != 1) {
    throw std::runtime_error(
        "MS has multiple fields, this does not work with the LOFAR beam "
        "library.");
  }

  if (common::HasColumn(fieldTable, "LOFAR_TILE_BEAM_DIR")) {
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
