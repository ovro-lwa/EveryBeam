#include "dish.h"
#include "../griddedresponse/dishgrid.h"

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

using everybeam::griddedresponse::DishGrid;
using everybeam::griddedresponse::GriddedResponse;
using everybeam::telescope::Dish;

Dish::Dish(casacore::MeasurementSet &ms, const Options &options)
    : Telescope(ms, options) {
  casacore::MSField field_table = ms.field();
  casacore::ArrayColumn<double> pointing_dir_col(
      field_table, casacore::MSField::columnName(casacore::MSField::DELAY_DIR));

  for (std::size_t field_id = 0; field_id != field_table.nrow(); ++field_id) {
    casacore::Array<double> pdir = pointing_dir_col(field_id);
    double pdir_ra = *pdir.cbegin();
    double pdir_dec = *(pdir.cbegin() + 1);
    ms_properties_.field_pointing.emplace_back(pdir_ra, pdir_dec);
  }
}

std::unique_ptr<GriddedResponse> Dish::GetGriddedResponse(
    const coords::CoordinateSystem &coordinate_system) {
  // Get and return GriddedResponse ptr
  std::unique_ptr<GriddedResponse> grid(new DishGrid(this, coordinate_system));
  return grid;
};