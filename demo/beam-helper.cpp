// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "beam-helper.h"

#include <limits>

#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Quanta/MVDirection.h>
#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/TableMeasures/ArrayQuantColumn.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/ms/MeasurementSets/MSAntenna.h>
#include <casacore/ms/MeasurementSets/MSAntennaEnums.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <aocommon/imagecoordinates.h>
#include <fitsio.h>

#include "./../cpp/common/mathutils.h"
#include "./../cpp/coords/itrfconverter.h"

using everybeam::vector2r_t;
using everybeam::vector3r_t;

void GetPhaseCentreInfo(casacore::MeasurementSet& ms, size_t fieldId,
                        double& ra, double& dec) {
  casacore::MSAntenna aTable = ms.antenna();
  size_t antennaCount = aTable.nrow();
  if (antennaCount == 0) throw std::runtime_error("No antennae in set");
  casacore::MPosition::ScalarColumn antPosColumn(
      aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
  casacore::MPosition ant1Pos = antPosColumn(0);
  casacore::MEpoch::ScalarColumn timeColumn(
      ms, ms.columnName(casacore::MSMainEnums::TIME));
  casacore::MSField fTable(ms.field());
  casacore::MDirection::ScalarColumn phaseDirColumn(
      fTable, fTable.columnName(casacore::MSFieldEnums::PHASE_DIR));
  casacore::MDirection phaseDir = phaseDirColumn(fieldId);
  casacore::MEpoch curtime = timeColumn(0);
  casacore::MeasFrame frame(ant1Pos, curtime);
  casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
  casacore::MDirection j2000 =
      casacore::MDirection::Convert(phaseDir, j2000Ref)();
  casacore::Vector<casacore::Double> j2000Val = j2000.getValue().get();
  ra = j2000Val[0];
  dec = j2000Val[1];
}

void GetThetaPhiDirectionsZenith(vector2r_t* thetaPhiDirections,
                                 size_t subgrid_size) {
  for (unsigned y = 0; y < subgrid_size; y++) {
    for (unsigned x = 0; x < subgrid_size; x++) {
      // Scale y,x to the interval -1, 1
      double l = (((double)2 * x) - subgrid_size) / subgrid_size;
      double m = (((double)2 * y) - subgrid_size) / subgrid_size;
      double n = sqrt(1 - l * l - m * m);

      // Compute direction in theta, phi
      vector2r_t direction_thetaphi;
      if (std::isfinite(n)) {
        // Convert direction to theta, phi
        vector3r_t direction_xyz = {(double)l, (double)m, n};
        direction_thetaphi = everybeam::cart2thetaphi(direction_xyz);
      } else {
        auto nan = std::numeric_limits<double>::quiet_NaN();
        direction_thetaphi = {nan, nan};
      }

      // Set direction
      thetaPhiDirections[y * subgrid_size + x] = direction_thetaphi;
    }
  }
}

void GetITRFDirections(vector3r_t* itrfDirections, size_t subgrid_size,
                       float image_size, double time, double ra, double dec) {
  float subgrid_pixelsize = image_size / subgrid_size;

  for (size_t y = 0; y < subgrid_size; y++) {
    for (size_t x = 0; x < subgrid_size; x++) {
      // IDG uses a flipped coordinate system which is moved by half a pixel:
      double dl = -subgrid_pixelsize;
      double dm = -subgrid_pixelsize;
      double pdl = -0.5 * dl;
      double pdm = 0.5 * dm;

      double l, m, n;

      aocommon::ImageCoordinates::XYToLM<double>(x, y, dl, dm, subgrid_size,
                                                 subgrid_size, l, m);

      l += pdl;
      m += pdm;
      n = sqrt(1.0 - l * l - m * m);

      const everybeam::coords::ItrfConverter itrf_converter(time);

      const vector3r_t l_vector_itrf =
          itrf_converter.RaDecToItrf(ra + M_PI / 2.0, 0);
      const vector3r_t m_vector_itrf =
          itrf_converter.RaDecToItrf(ra, dec + M_PI / 2.0);
      const vector3r_t n_vector_itrf = itrf_converter.RaDecToItrf(ra, dec);

      const vector3r_t itrf_direction{
          l * l_vector_itrf[0] + m * m_vector_itrf[0] + n * n_vector_itrf[0],
          l * l_vector_itrf[1] + m * m_vector_itrf[1] + n * n_vector_itrf[1],
          l * l_vector_itrf[2] + m * m_vector_itrf[2] + n * n_vector_itrf[2]};

      itrfDirections[y * subgrid_size + x] = itrf_direction;
    }
  }
}

void CheckFitsStatus(int status, const std::string& filename) {
  if (status) {
    /* fits_get_errstatus returns at most 30 characters */
    char err_text[31];
    fits_get_errstatus(status, err_text);
    char err_msg[81];
    std::stringstream errMsg;
    errMsg << "CFITSIO reported error when performing IO on file '" << filename
           << "': " << err_text << " (";
    while (fits_read_errmsg(err_msg)) {
      errMsg << err_msg;
    }
    errMsg << ')';
    throw std::runtime_error(errMsg.str());
  }
}

template <typename NumType>
void WriteFits(const std::string& filename, const NumType* image, size_t width,
               size_t height) {
  // Open file
  fitsfile* fptr;
  int status = 0;
  fits_create_file(&fptr, (std::string("!") + filename).c_str(), &status);
  CheckFitsStatus(status, filename);

  // Write header
  long axes[] = {(long)width, (long)height};
  fits_create_img(fptr, FLOAT_IMG, 2, axes, &status);
  CheckFitsStatus(status, filename);

  // Write image
  long pixel[4] = {1, 1, 1, 1};
  double nullValue = std::numeric_limits<double>::max();
  size_t totalSize = width * height;
  fits_write_pixnull(fptr, TDOUBLE, pixel, totalSize, (void*)image, &nullValue,
                     &status);
  CheckFitsStatus(status, filename);

  // Close file
  fits_close_file(fptr, &status);
  CheckFitsStatus(status, filename);
}

void StoreBeam(const std::string& filename, const std::complex<float>* buffer,
               size_t nStations, size_t width, size_t height) {
  size_t ny = floor(sqrt(nStations)), nx = (nStations + ny - 1) / ny;
  std::cout << "Storing " << filename << " (" << nStations << " ant, " << nx
            << " x " << ny << ")\n";
  std::vector<double> img(nx * ny * width * height, 0.0);
  for (size_t ant = 0; ant != nStations; ++ant) {
    size_t xCorner = (ant % nx) * width, yCorner = (ant / nx) * height;
    for (size_t y = 0; y != height; ++y) {
      for (size_t x = 0; x != width; ++x) {
        std::complex<float> response = 0;
        for (auto pol = 0; pol < 4; pol++) {
          response += *buffer * std::conj(*buffer);
          ++buffer;
        }
        img[(yCorner + y) * width * nx + x + xCorner] = abs(response) / 2;
      }
    }
  }
  WriteFits<double>(filename, img.data(), nx * width, ny * height);
}

void GetRaDecZenith(vector3r_t position, double time, double& ra, double& dec) {
  casacore::MEpoch timeEpoch(casacore::Quantity(time, "s"));
  casacore::MVPosition mvPosition(position[0], position[1], position[2]);
  casacore::MPosition mPosition(mvPosition, casacore::MPosition::ITRF);
  casacore::MeasFrame mFrame(timeEpoch, mPosition);
  auto converter = casacore::MDirection::Convert(
      casacore::MDirection::ITRF,
      casacore::MDirection::Ref(casacore::MDirection::J2000, mFrame));
  casacore::MVDirection mvZenith(position[0], position[1], position[2]);
  casacore::MDirection mZenith(mvZenith, casacore::MDirection::ITRF);
  casacore::MVDirection zenithRaDec = converter(mZenith).getValue();
  ra = zenithRaDec.getAngle().getValue()[0];
  dec = zenithRaDec.getAngle().getValue()[1];
}

std::string GetFieldName(casacore::MeasurementSet& ms, unsigned int field_id) {
  casacore::Table fieldTable = ms.keywordSet().asTable("LOFAR_ANTENNA_FIELD");
  casacore::ROScalarColumn<casacore::String> nameColumn(fieldTable, "NAME");
  return nameColumn(field_id);
}

std::string GetStationName(casacore::MeasurementSet& ms,
                           unsigned int station_id) {
  casacore::Table antennaTable = ms.keywordSet().asTable("ANTENNA");
  casacore::ROScalarColumn<casacore::String> nameColumn(antennaTable, "NAME");
  return nameColumn(station_id);
}

unsigned int GetNrAntennas(casacore::MeasurementSet& ms,
                           unsigned int field_id) {
  casacore::Table fieldTable = ms.keywordSet().asTable("LOFAR_ANTENNA_FIELD");
  casacore::ArrayQuantColumn<double> offsetColumn(fieldTable, "ELEMENT_OFFSET",
                                                  "m");
  casacore::Matrix<casacore::Quantity> aips_offset = offsetColumn(field_id);
  return aips_offset.ncolumn();
}
