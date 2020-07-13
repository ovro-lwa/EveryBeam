#include "lofar.h"
#include "./../gridded_response/lofargrid.h"
#include "./../common/math_utils.h"
#include "./../common/casa_utils.h"

#include <aocommon/banddata.h>
#include <cassert>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using namespace everybeam;
using namespace everybeam::telescope;
using namespace casacore;

// LOFAR Telescope specific utils
namespace {

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

// TODO: move to common utils?
vector3r_t transformToFieldCoordinates(
    const vector3r_t &position, const Antenna::CoordinateSystem::Axes &axes) {
  const vector3r_t result{dot(position, axes.p), dot(position, axes.q),
                          dot(position, axes.r)};
  return result;
}

// // Seems LOFAR specific?
// void transformToFieldCoordinates(TileConfig &config,
//                                  const Antenna::CoordinateSystem::Axes &axes)
//                                  {
//   for (unsigned int i = 0; i < config.size(); ++i) {
//     const vector3r_t position = config[i];
//     config[i][0] = dot(position, axes.p);
//     config[i][1] = dot(position, axes.q);
//     config[i][2] = dot(position, axes.r);
//   }
// }

// LOFAR specific? Or is this generic for each telescope?
BeamFormer::Ptr readAntennaField(const Table &table, std::size_t id,
                                 ElementResponse::Ptr element_response) {
  Antenna::CoordinateSystem coordinate_system =
      common::readCoordinateSystem(table, id);

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

  //     TileConfig tile_config;
  //     if(name != "LBA") readTileConfig(table, id);
  //     transformToFieldCoordinates(tile_config, coordinate_system.axes);

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
      antenna = Element::Ptr(
          new Element(antenna_coordinate_system, element_response, id));

      // TODO
      // HBA, HBA0, HBA1
      // antenna = make_tile(name, coordinate_system, tile_config,
      // element_response);
    }

    antenna->m_enabled[0] = !aips_flag(0, i);
    antenna->m_enabled[1] = !aips_flag(1, i);
    beam_former->AddAntenna(antenna);
  }
  return beam_former;
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
}  // namespace

LOFAR::LOFAR(MeasurementSet &ms, const ElementResponseModel model,
             const Options &options)
    : Telescope(ms, model, options) {
  ReadAllStations(ms, model);

  // Populate MeasurementSet properties struct
  aocommon::BandData band(ms.spectralWindow());
  MSAntenna antenna(ms.antenna());
  MPosition::ScalarColumn antenna_pos_col(
      antenna, antenna.columnName(MSAntennaEnums::POSITION));
  MEpoch::ScalarColumn time_column(ms,
                                   ms.columnName(casacore::MSMainEnums::TIME));

  // Following is ms.field() related, first check whether field complies with
  // LOFAR field
  if (ms.field().nrow() != 1)
    throw std::runtime_error("LOFAR MeasurementSet has multiple fields");

  if (!ms.field().tableDesc().isColumn("LOFAR_TILE_BEAM_DIR")) {
    throw std::runtime_error("LOFAR_TILE_BEAM_DIR column not found");
  }

  casacore::ScalarMeasColumn<casacore::MDirection> delay_dir_col(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));

  casacore::ArrayMeasColumn<casacore::MDirection> tile_beam_dir_col(
      ms.field(), "LOFAR_TILE_BEAM_DIR");

  // Populate struct
  ms_properties_ = {.subband_freq = band.CentreFrequency(),
                    .delay_dir = delay_dir_col(0),
                    .tile_beam_dir = *(tile_beam_dir_col(0).data())};
}

std::unique_ptr<gridded_response::GriddedResponse> LOFAR::GetGriddedResponse(
    const coords::CoordinateSystem &coordinate_system) {
  // Get and return GriddedResponse ptr
  std::unique_ptr<gridded_response::GriddedResponse> grid(
      new gridded_response::LOFARGrid(this, coordinate_system));
  // gridded_response::GriddedResponse grid(LOFARGrid(this, coordinate_system));
  return grid;
};

Station::Ptr LOFAR::ReadStation(const MeasurementSet &ms, std::size_t id,
                                const ElementResponseModel model) const {
  ROMSAntennaColumns antenna(ms.antenna());
  assert(nstations_ > id && !antenna.flagRow()(id));

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

  Table tab_field = common::getSubTable(ms, "LOFAR_ANTENNA_FIELD");
  tab_field = tab_field(tab_field.col("ANTENNA_ID") == static_cast<Int>(id));

  // The Station will consist of a BeamFormer that combines the fields
  // coordinate system is ITRF
  // phase reference is station position
  auto beam_former = BeamFormer::Ptr(new BeamFormer(
      Antenna::IdentityCoordinateSystem, station->phaseReference()));

  for (std::size_t i = 0; i < tab_field.nrow(); ++i) {
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

  // TODO: rotate coordinate system for antenna
  auto element = Element::Ptr(new Element(
      coordinate_system, station->get_element_response(), element_id));
  station->SetElement(element);

  return station;
}