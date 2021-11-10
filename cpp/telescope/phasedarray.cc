// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "phasedarray.h"
#include "../common/casautils.h"

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

namespace everybeam {
namespace telescope {
void PhasedArray::CalculatePreappliedBeamOptions(
    const casacore::MeasurementSet &ms, const std::string &data_column_name,
    casacore::MDirection &preapplied_beam_dir,
    CorrectionMode &correction_mode) {
  casacore::ScalarMeasColumn<casacore::MDirection> referenceDirColumn(
      ms.field(),
      casacore::MSField::columnName(casacore::MSFieldEnums::REFERENCE_DIR));
  preapplied_beam_dir = referenceDirColumn(0);

  // Read beam keywords of input datacolumn
  // NOTE: it is somewhat confusing that "LOFAR" is explicitly mentioned
  // in the keyword names. The keywords, however, naturally extend to the
  // OSKAR telescope as well.
  casacore::ArrayColumn<std::complex<float>> dataCol(ms, data_column_name);
  if (dataCol.keywordSet().isDefined("LOFAR_APPLIED_BEAM_MODE")) {
    correction_mode = ParseCorrectionMode(
        dataCol.keywordSet().asString("LOFAR_APPLIED_BEAM_MODE"));
    switch (correction_mode) {
      case CorrectionMode::kNone:
        break;
      case CorrectionMode::kElement:
      case CorrectionMode::kArrayFactor:
      case CorrectionMode::kFull:
        casacore::String error;
        casacore::MeasureHolder mHolder;
        if (!mHolder.fromRecord(error, dataCol.keywordSet().asRecord(
                                           "LOFAR_APPLIED_BEAM_DIR"))) {
          throw std::runtime_error(
              "Error while reading LOFAR_APPLIED_BEAM_DIR keyword: " + error);
        }
        preapplied_beam_dir = mHolder.asMDirection();
        break;
    }
  } else {
    correction_mode = CorrectionMode::kNone;
  }
}
}  // namespace telescope
}  // namespace everybeam