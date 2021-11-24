// msreadutils.cc: Utility functions to read the meta data relevant for
// simulating the beam from phased array (LOFAR/OSKAR) observations stored in MS
// format.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "msreadutils.h"
#include "load.h"
#include "beamformeridenticalantennas.h"
#include "beamformerlofarhba.h"
#include "beamformerlofarlba.h"
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

using casacore::ArrayMeasColumn;
using casacore::ArrayQuantColumn;
using casacore::Bool;
using casacore::Double;
using casacore::Int;
using casacore::Matrix;
using casacore::MDirection;
using casacore::MeasurementSet;
using casacore::MPosition;
using casacore::MVPosition;
using casacore::Quantity;
using casacore::ROArrayColumn;
using casacore::ROMSAntennaColumns;
using casacore::ROScalarColumn;
using casacore::ROScalarMeasColumn;
using casacore::String;
using casacore::Table;

namespace everybeam {
namespace {
using TileConfig = std::array<vector3r_t, 16>;

constexpr Antenna::CoordinateSystem::Axes oskar_antenna_orientation = {
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

constexpr Antenna::CoordinateSystem::Axes lofar_antenna_orientation = {
    {
        -M_SQRT1_2,
        -M_SQRT1_2,
        0.0,
    },
    {
        M_SQRT1_2,
        -M_SQRT1_2,
        0.0,
    },
    {0.0, 0.0, 1.0},
};

TileConfig ReadTileConfig(const Table &table, unsigned int row) {
  ArrayQuantColumn<Double> c_tile_offset(table, "TILE_ELEMENT_OFFSET", "m");

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
  for (auto &val : config) {
    const vector3r_t position = val;
    val[0] = dot(position, axes.p);
    val[1] = dot(position, axes.q);
    val[2] = dot(position, axes.r);
  }
}

vector3r_t TransformToFieldCoordinates(
    const vector3r_t &position, const Antenna::CoordinateSystem::Axes &axes) {
  const vector3r_t result{dot(position, axes.p), dot(position, axes.q),
                          dot(position, axes.r)};
  return result;
}

std::shared_ptr<BeamFormer> MakeTile(const vector3r_t &position,
                                     const TileConfig &tile_config,
                                     ElementResponse::Ptr element_response) {
  std::shared_ptr<BeamFormer> tile =
      std::make_shared<BeamFormerIdenticalAntennas>(position);

  for (unsigned int id = 0; id < tile_config.size(); id++) {
    vector3r_t antenna_position = tile_config[id];

    Antenna::CoordinateSystem antenna_coordinate_system;
    antenna_coordinate_system.origin = antenna_position;
    antenna_coordinate_system.axes = lofar_antenna_orientation;

    std::shared_ptr<Antenna> antenna = std::make_shared<Element>(
        antenna_coordinate_system, element_response, id);
    tile->AddAntenna(antenna);
  }
  return tile;
}

// Make a dedicated HBA "Hamaker" tile, saving only one element, and 16
// element positions
void MakeTile(std::shared_ptr<BeamFormerLofarHBA> beamformer,
              const TileConfig &tile_config,
              ElementResponse::Ptr element_response) {
  for (unsigned int id = 0; id < tile_config.size(); id++) {
    vector3r_t antenna_position = tile_config[id];

    Antenna::CoordinateSystem antenna_coordinate_system;
    antenna_coordinate_system.origin = antenna_position;
    antenna_coordinate_system.axes = lofar_antenna_orientation;

    std::shared_ptr<ElementHamaker> antenna = std::make_shared<ElementHamaker>(
        antenna_coordinate_system, element_response, id);

    // Only element 1 needs to be stored as an element
    if (id == 0) {
      beamformer->SetElement(antenna);
    }
    // All positions need to be stored, however
    beamformer->AddElementPosition(antenna_position);
  }
}

std::shared_ptr<Antenna> ReadAntennaFieldLofar(
    const Table &table, unsigned int id,
    ElementResponse::Ptr element_response) {
  Antenna::CoordinateSystem coordinate_system =
      common::ReadCoordinateSystem(table, id);

  ScalarColumn<String> c_name(table, "NAME");
  ArrayQuantColumn<Double> c_offset(table, "ELEMENT_OFFSET", "m");
  ArrayColumn<Bool> c_flag(table, "ELEMENT_FLAG");

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
  if (element_response->GetModel() == ElementResponseModel::kHamaker) {
    if (name != "LBA") {
      // Then HBA, HBA0 or HBA1
      beam_former = std::make_shared<BeamFormerLofarHBA>(coordinate_system);
    } else {
      beam_former = std::make_shared<BeamFormerLofarLBA>(coordinate_system);
    }
  } else if (element_response->GetModel() == ElementResponseModel::kLOBES) {
    // BeamFormer is assigned a FieldResponse, for which common field quantities
    // can be precomputed
    beam_former = std::make_shared<BeamFormer>(
        coordinate_system,
        std::dynamic_pointer_cast<FieldResponse>(*element_response.get()));
  } else {
    // All Tiles / elements should be kept unique, so work with generic
    // BeamFormer
    beam_former = std::make_shared<BeamFormer>(coordinate_system);
  }

  TransformToFieldCoordinates(tile_config, coordinate_system.axes);

  for (size_t i = 0; i < aips_offset.ncolumn(); ++i) {
    vector3r_t antenna_position = {aips_offset(0, i).getValue(),
                                   aips_offset(1, i).getValue(),
                                   aips_offset(2, i).getValue()};
    antenna_position =
        TransformToFieldCoordinates(antenna_position, coordinate_system.axes);
    std::shared_ptr<Antenna> antenna;
    Antenna::CoordinateSystem antenna_coordinate_system{
        antenna_position, lofar_antenna_orientation};

    if (name == "LBA") {
      antenna = std::make_shared<Element>(antenna_coordinate_system,
                                          element_response, i);
      if (element_response->GetModel() == kHamaker) {
        // Cast to LOFAR LBA
        std::shared_ptr<BeamFormerLofarLBA> beam_former_lba =
            std::static_pointer_cast<BeamFormerLofarLBA>(beam_former);
        // Store only one Element
        if (i == 0) {
          std::shared_ptr<Element> element =
              std::dynamic_pointer_cast<Element>(antenna);
          beam_former_lba->SetElement(element);
        }
        // Store Element position and Element enabled in x/y?
        beam_former_lba->AddElementPosition(antenna_position);
        beam_former_lba->AddElementEnabled(
            std::array<bool, 2>{!aips_flag(0, i), !aips_flag(0, i)});
      } else {
        antenna->enabled_[0] = !aips_flag(0, i);
        antenna->enabled_[1] = !aips_flag(1, i);

        std::shared_ptr<BeamFormer> beam_former_lba =
            std::static_pointer_cast<BeamFormer>(beam_former);
        beam_former_lba->AddAntenna(antenna);
      }
    } else {
      // name is HBA, HBA0 or HBA1
      if (element_response->GetModel() == kHamaker) {
        std::shared_ptr<BeamFormerLofarHBA> beam_former_hba =
            std::static_pointer_cast<BeamFormerLofarHBA>(beam_former);

        // Tile positions are unique
        beam_former_hba->AddTilePosition(antenna_position);

        // Store only one tile
        if (i == 0) {
          MakeTile(beam_former_hba, tile_config, element_response);
        }
        // Tile enabled in x/y?
        beam_former_hba->AddTileEnabled(
            std::array<bool, 2>{!aips_flag(0, i), !aips_flag(0, i)});
      } else {
        std::shared_ptr<BeamFormerIdenticalAntennas> beam_former_hba =
            std::static_pointer_cast<BeamFormerIdenticalAntennas>(beam_former);
        antenna = MakeTile(antenna_position, tile_config, element_response);
        antenna->enabled_[0] = !aips_flag(0, i);
        antenna->enabled_[1] = !aips_flag(1, i);
        beam_former_hba->AddAntenna(antenna);
      }
    }
  }
  return beam_former;
}

std::shared_ptr<Element> AartfaacElement(
    const MeasurementSet &ms, size_t station_id,
    ElementResponse::Ptr element_response) {
  Table table = common::GetSubTable(ms, "ANTENNA");

  ScalarColumn<String> antenna_type_column(ms.observation(),
                                           everybeam::kAartfaacAntennaTypeName);
  const std::string ant_type = antenna_type_column(0);

  Antenna::CoordinateSystem coordinate_system =
      common::ReadAartfaacCoordinateSystem(table, station_id);

  const size_t id = 0;
  std::shared_ptr<Element> antenna =
      std::make_shared<Element>(coordinate_system, element_response, id);
  antenna->enabled_[0] = true;
  antenna->enabled_[1] = true;
  return antenna;
}

std::shared_ptr<BeamFormer> ReadAntennaFieldMSv3(
    const Table &table, size_t station_id,
    ElementResponse::Ptr element_response) {
  Antenna::CoordinateSystem coordinate_system =
      common::ReadCoordinateSystem(table, station_id);
  std::shared_ptr<BeamFormer> beam_former =
      std::make_shared<BeamFormerIdenticalAntennas>(coordinate_system);

  ArrayQuantColumn<Double> c_offset(table, "ELEMENT_OFFSET", "m");
  ArrayColumn<Bool> c_flag(table, "ELEMENT_FLAG");

  // Read element offsets and flags.
  Matrix<Quantity> aips_offset = c_offset(station_id);

  assert(aips_offset.shape().isEqual(IPosition(2, aips_offset.nrow(), 3)));

  Matrix<Bool> aips_flag = c_flag(station_id);
  assert(aips_flag.shape().isEqual(IPosition(2, aips_offset.nrow(), 2)));

  for (size_t i = 0; i < aips_offset.nrow(); ++i) {
    vector3r_t antenna_position = {aips_offset(i, 0).getValue(),
                                   aips_offset(i, 1).getValue(),
                                   aips_offset(i, 2).getValue()};
    antenna_position =
        TransformToFieldCoordinates(antenna_position, coordinate_system.axes);
    std::shared_ptr<Antenna> antenna;
    Antenna::CoordinateSystem antenna_coordinate_system{
        antenna_position, oskar_antenna_orientation};
    antenna = std::make_shared<Element>(antenna_coordinate_system,
                                        element_response, i);

    antenna->enabled_[0] = !aips_flag(i, 0);
    antenna->enabled_[1] = !aips_flag(i, 1);
    beam_former->AddAntenna(antenna);
  }
  return beam_former;
}

std::shared_ptr<BeamFormer> LofarStationBeamFormer(
    const MeasurementSet &ms, size_t station_id,
    const vector3r_t &phase_reference, ElementResponse::Ptr element_response) {
  std::shared_ptr<BeamFormer> beam_former;

  Table tab_field = common::GetSubTable(ms, "LOFAR_ANTENNA_FIELD");
  tab_field =
      tab_field(tab_field.col("ANTENNA_ID") == static_cast<Int>(station_id));

  // The Station will consist of a BeamFormer that combines the fields
  // coordinate system is ITRF
  // phase reference is station position
  beam_former = std::make_shared<BeamFormer>(Antenna::IdentityCoordinateSystem,
                                             phase_reference);

  for (size_t i = 0; i < tab_field.nrow(); ++i) {
    beam_former->AddAntenna(
        ReadAntennaFieldLofar(tab_field, i, element_response));
  }

  // TODO
  // If There is only one field, the top level beamformer is not needed
  // and the station antenna can be set to the beamformer of the field
  // station->SetAntenna(beam_former);
  return beam_former;
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

std::shared_ptr<BeamFormer> MSv3StationBeamFormer(
    const MeasurementSet &ms, size_t station_id,
    ElementResponse::Ptr element_response) {
  Table tab_phased_array = common::GetSubTable(ms, "PHASED_ARRAY");

  // The Station will consist of a BeamFormer that combines the fields
  // coordinate system is ITRF
  auto beam_former =
      ReadAntennaFieldMSv3(tab_phased_array, station_id, element_response);
  return beam_former;
}
}  // namespace

std::shared_ptr<Station> ReadSingleStation(const casacore::MeasurementSet &ms,
                                           unsigned int id,
                                           const Options &options) {
  TelescopeType telescope_type = GetTelescopeType(ms);
  if (telescope_type != TelescopeType::kLofarTelescope &&
      telescope_type != TelescopeType::kAARTFAAC &&
      telescope_type != TelescopeType::kOSKARTelescope) {
    throw std::runtime_error(
        "MSReadUtils found an unknown telescope type in MS, return value of "
        "GetTelescopeType(ms) should be one of "
        "kLofarTelescope, kAARTFAAC, or kOSKARTelescope.");
  }

  ROMSAntennaColumns antenna(ms.antenna());
  assert(antenna.nrow() > id && !antenna.flagRow()(id));

  // Get station name.
  const std::string name(antenna.name()(id));

  // Get station position (ITRF).
  MPosition mPosition =
      MPosition::Convert(antenna.positionMeas()(id), MPosition::ITRF)();
  MVPosition mvPosition = mPosition.getValue();
  const vector3r_t position = {mvPosition(0), mvPosition(1), mvPosition(2)};

  // Create station
  std::shared_ptr<Station> station =
      std::make_shared<Station>(name, position, options);

  // Set the top level beamformer (that might contain nested beam formers)
  if (telescope_type == TelescopeType::kOSKARTelescope) {
    // OSKAR telescope
    auto beam_former =
        MSv3StationBeamFormer(ms, id, station->GetElementResponse());
    station->SetAntenna(beam_former);
  } else {
    // LOFAR Telescope or AARTFAAC
    station->SetPhaseReference(ReadStationPhaseReference(ms.antenna(), id));
    if (telescope_type == TelescopeType::kLofarTelescope) {
      auto beam_former = LofarStationBeamFormer(
          ms, id, station->GetPhaseReference(), station->GetElementResponse());
      station->SetAntenna(beam_former);
    } else if (telescope_type == TelescopeType::kAARTFAAC) {
      auto element = AartfaacElement(ms, id, station->GetElementResponse());
      station->SetAntenna(element);
    }
  }
  return station;
}

MDirection ReadTileBeamDirection(const casacore::MeasurementSet &ms) {
  TelescopeType telescope_type = GetTelescopeType(ms);
  if (telescope_type != TelescopeType::kLofarTelescope &&
      telescope_type != TelescopeType::kAARTFAAC) {
    throw std::runtime_error(
        "Tile beam direction requested. This does not work with MS other than "
        "LOFAR or AARTFAAC.");
  }

  MDirection tile_beam_dir;

  Table field_table = common::GetSubTable(ms, "FIELD");

  if (field_table.nrow() != 1) {
    throw std::runtime_error(
        "MS has multiple fields, this does not work with the LOFAR beam "
        "library.");
  }

  if (common::HasColumn(field_table, "LOFAR_TILE_BEAM_DIR")) {
    ROArrayMeasColumn<MDirection> tile_beam_col(field_table,
                                                "LOFAR_TILE_BEAM_DIR");
    tile_beam_dir = *(tile_beam_col(0).data());
  } else {
    ROArrayMeasColumn<MDirection> tile_beam_col(field_table, "DELAY_DIR");
    tile_beam_dir = *(tile_beam_col(0).data());
  }

  return tile_beam_dir;
}

}  // namespace everybeam