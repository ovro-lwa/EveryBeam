// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <complex>
#include <cstddef>
#include <string>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "../cpp/common/types.h"

[[gnu::visibility("default")]] void GetPhaseCentreInfo(
    casacore::MeasurementSet& ms, size_t fieldId, double& ra, double& dec);

[[gnu::visibility("default")]] void GetThetaPhiDirectionsZenith(
    everybeam::vector2r_t* raDecDirections, size_t subgrid_size);

[[gnu::visibility("default")]] void GetITRFDirections(
    everybeam::vector3r_t* itrfDirections, size_t subgrid_size,
    float image_size, double time, double phaseCentreRA, double phaseCentreDec);

[[gnu::visibility("default")]] void StoreBeam(const std::string& filename,
                                              const std::complex<float>* buffer,
                                              size_t nStations, size_t width,
                                              size_t height);

[[gnu::visibility("default")]] void GetRaDecZenith(
    everybeam::vector3r_t position, double time, double& ra, double& dec);

[[gnu::visibility("default")]] std::string GetFieldName(
    casacore::MeasurementSet& ms, unsigned int field_id = 0);

[[gnu::visibility("default")]] std::string GetStationName(
    casacore::MeasurementSet& ms, unsigned int station_id);

[[gnu::visibility("default")]] unsigned int GetNrAntennas(
    casacore::MeasurementSet& ms, unsigned int field_id);