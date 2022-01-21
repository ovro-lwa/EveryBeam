// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>

#include "config.h"
#include "lobeselementresponse.h"

#include "../common/sphericalharmonics.h"

#include <map>
#include <tuple>
#include <complex>
#include <H5Cpp.h>
#include <iostream>
#include <boost/filesystem.hpp>
#include <string>

namespace everybeam {

namespace {
/**
 * @brief Search for LOBES h5 coefficient file
 * on the suggested path \param search_path. Returns
 * an empty string if the file cannot be found.
 *
 * @param search_path Search path
 * @param station_name Station name, as read from MS
 * @return std::string Path to file or empty string if file cannot be found
 */

boost::filesystem::path FindCoeffFile(const std::string &search_path,
                                      const std::string &station_name) {
  std::string station_file = "LOBES_" + station_name + ".h5";

  boost::filesystem::path p =
      search_path.empty()
          ? boost::filesystem::path(std::string{EVERYBEAM_DATA_DIR} +
                                    std::string{"/lobes"}) /
                station_file
          : boost::filesystem::path(search_path) / station_file;
  return p;
}
}  // namespace

LOBESElementResponse::LOBESElementResponse(const std::string &name,
                                           const Options &options) {
  boost::filesystem::path coeff_file_path =
      FindCoeffFile(options.coeff_path, name);
  H5::H5File h5file;

  if (!boost::filesystem::exists(coeff_file_path)) {
    throw std::runtime_error("LOBES coeffcients file: " +
                             coeff_file_path.string() + " does not exists");
  }

  try {
    h5file.openFile(coeff_file_path.c_str(), H5F_ACC_RDONLY);
  } catch (const H5::FileIException &e) {
    throw std::runtime_error("Could not open LOBES coeffcients file: " +
                             coeff_file_path.string());
  }

  H5::DataSet dataset = h5file.openDataSet("coefficients");
  H5::DataSpace dataspace = dataset.getSpace();
  int nr_elements = dataspace.getSimpleExtentNpoints();

  const std::string REAL("r");
  const std::string IMAG("i");

  // Create HDF5 complex datatype
  H5::CompType h5_dcomplex(sizeof(std::complex<double>));
  h5_dcomplex.insertMember(REAL, 0, H5::PredType::NATIVE_DOUBLE);
  h5_dcomplex.insertMember(IMAG, sizeof(double), H5::PredType::NATIVE_DOUBLE);

  // Get the number of dimensions in the dataspace.
  int ndims_coefficients = dataspace.getSimpleExtentNdims();

  // Get the dimension size of each dimension in the dataspace and display them.
  std::vector<hsize_t> dims_coefficients(ndims_coefficients);
  dataspace.getSimpleExtentDims(dims_coefficients.data(), nullptr);
  coefficients_shape_ = std::vector<unsigned int>(dims_coefficients.begin(),
                                                  dims_coefficients.end());

  // Read coefficients
  coefficients_.resize(coefficients_shape_[0], coefficients_shape_[1],
                       coefficients_shape_[2], coefficients_shape_[3]);
  dataset.read(coefficients_.data(), h5_dcomplex);

  // Frequencies
  dataset = h5file.openDataSet("frequencies");
  dataspace = dataset.getSpace();
  nr_elements = dataspace.getSimpleExtentNpoints();

  frequencies_.resize(nr_elements);
  dataset.read(frequencies_.data(), H5::PredType::NATIVE_DOUBLE);

  // nms
  dataset = h5file.openDataSet("nms");
  dataspace = dataset.getSpace();
  nr_elements = dataspace.getSimpleExtentNpoints();

  // Get the number of dimensions in the dataspace.
  int ndims_nms = dataspace.getSimpleExtentNdims();

  // Get the dimension size of each dimension in the dataspace and display them.
  std::vector<hsize_t> dims_nms(ndims_nms);
  dataspace.getSimpleExtentDims(dims_nms.data(), nullptr);

  nms_.resize(dims_nms[0]);
  dataset.read(nms_.data(), H5::PredType::NATIVE_INT);
}

LOBESElementResponse::BaseFunctions LOBESElementResponse::ComputeBaseFunctions(
    double theta, double phi) const {
  LOBESElementResponse::BaseFunctions base_functions(nms_.size(), 2);
  base_functions.setZero();

  for (size_t i = 0; i < nms_.size(); ++i) {
    auto nms = nms_[i];
    std::complex<double> q2, q3;
    std::tie(q2, q3) =
        everybeam::common::F4far_new(nms.s, nms.m, nms.n, theta, phi);
    base_functions(i, 0) = q2;
    base_functions(i, 1) = q3;
  }
  return base_functions;
}

aocommon::MC2x2 LOBESElementResponse::Response(int element_id, double freq,
                                               double theta, double phi) const {
  // Clip directions below the horizon.
  if (theta >= M_PI_2) {
    return aocommon::MC2x2::Zero();
  }

  bool clear_basefunctions = false;
  if (basefunctions_.rows() == 0) {
    // Fill basefunctions if not cached at this stage (via SetFieldQuantities).
    // Also, set clear_basefunctions to true, to clear the basefunctions_ at the
    // end of this response calculation.
    basefunctions_ = ComputeBaseFunctions(theta, phi);
    clear_basefunctions = true;
  }

  const int freq_idx = FindFrequencyIdx(freq);
  std::complex<double> xx = {0}, xy = {0}, yx = {0}, yy = {0};

  const int nr_rows = basefunctions_.rows();
  if (nr_rows == 0) {
    throw std::runtime_error(
        "Number of rows in basefunctions_ member is 0. Did you run "
        "SetFieldQuantities?");
  }

  for (int i = 0; i < nr_rows; ++i) {
    const std::complex<double> q2 = basefunctions_(i, 0);
    const std::complex<double> q3 = basefunctions_(i, 1);
    xx += q2 * coefficients_(0, freq_idx, element_id, i);
    xy += q3 * coefficients_(0, freq_idx, element_id, i);
    yx += q2 * coefficients_(1, freq_idx, element_id, i);
    yy += q3 * coefficients_(1, freq_idx, element_id, i);
  }

  if (clear_basefunctions) {
    // Do a destructive resize
    basefunctions_.resize(0, 2);
  }
  return aocommon::MC2x2(xx, xy, yx, yy);
}

std::shared_ptr<LOBESElementResponse> LOBESElementResponse::GetInstance(
    const std::string &name, const Options &options) {
  static std::map<std::string, std::shared_ptr<LOBESElementResponse>>
      name_response_map;

  auto entry = name_response_map.find(name);
  if (entry == name_response_map.end()) {
    entry = name_response_map.insert(
        entry, {name, std::make_shared<LOBESElementResponse>(name, options)});
  }
  return entry->second;
}

}  // namespace everybeam
