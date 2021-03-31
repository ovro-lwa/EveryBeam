// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../cpp/elementresponse.h"
#include "../cpp/station.h"
#include "../cpp/msreadutils.h"
#include "../cpp/coords/coordutils.h"

#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/ms/MeasurementSets/MSAntenna.h>
#include <casacore/tables/Tables/TableKeyword.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>

using namespace everybeam;

void GetPhaseCentreInfo(casacore::MeasurementSet& ms, size_t fieldId,
                        double& ra, double& dec);

void GetThetaPhiDirectionsZenith(vector2r_t* raDecDirections,
                                 size_t subgrid_size);

void GetITRFDirections(vector3r_t* itrfDirections, size_t subgrid_size,
                       float image_size, double time, double phaseCentreRA,
                       double phaseCentreDec);

void StoreBeam(const std::string& filename, const std::complex<float>* buffer,
               size_t nStations, size_t width, size_t height);

void GetRaDecZenith(vector3r_t position, double time, double& ra, double& dec);

std::string GetFieldName(casacore::MeasurementSet& ms,
                         unsigned int field_id = 0);

std::string GetStationName(casacore::MeasurementSet& ms,
                           unsigned int station_id);

unsigned int GetNrAntennas(casacore::MeasurementSet& ms, unsigned int field_id);