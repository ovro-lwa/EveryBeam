// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <charconv>
#include <cmath>
#include <complex>
#include <filesystem>
#include <map>
#include <optional>
#include <string_view>

#include <aocommon/throwruntimeerror.h>
#include <boost/algorithm/string/predicate.hpp>
#include <H5Cpp.h>

#include "../common/mathutils.h"
#include "../common/sphericalharmonics.h"

#include "config.h"
#include "lobeselementresponse.h"
#include "lobeselementresponsefixeddirection.h"

// There are two main modi for the AARTFAAC telescope, AARTFAAC-6 and
// AARTFAAC-12. To properly use AARTFAAC in LOBEs mode the coefficients of all
// stations need to be available. At the moment of writing only a partial set
// is available. This means only AARTFAAC-6 is tested.
static const std::array<std::string_view, 12> kAartfaacStationNames{
    // Available
    "CS002LBA", "CS003LBA", "CS004LBA", "CS005LBA", "CS006LBA", "CS007LBA",
    "CS001LBA", "CS011LBA", "CS013LBA",
    // Currently unavailable
    "CS017LBA", "CS021LBA", "CS032LBA"};

struct AartfaacStation {
  std::string_view station;
  int element;
};

template <class T>
static T ExtractIntegral(std::string_view string) {
  int value;
  std::from_chars_result result =
      std::from_chars(string.begin(), string.end(), value);
  if (result.ec != std::errc{} || result.ptr != string.end()) {
    aocommon::ThrowRuntimeError("The value '", string,
                                "' can't be converted to a number");
  }
  return value;
}

enum class AartfaacElements { kInner, kOuter };

static std::optional<AartfaacStation> GetAartfaacStation(
    std::string_view station_name, AartfaacElements elements) {
  if (!boost::starts_with(station_name, "A12_")) {
    return {};
  }

  station_name.remove_prefix(4);
  const int id = ExtractIntegral<int>(station_name);
  const size_t station_id = id / 48;
  const int element_id =
      id % 48 + (elements == AartfaacElements::kInner ? 0 : 48);

  if (station_id >= kAartfaacStationNames.size()) {
    aocommon::ThrowRuntimeError("Aartfaac station id '", station_id,
                                "' is invalid");
  }
  return AartfaacStation{kAartfaacStationNames[station_id], element_id};
}

namespace everybeam {

static const H5::CompType kH5Dcomplex = [] {
  const std::string REAL("r");
  const std::string IMAG("i");

  H5::CompType h5_dcomplex(sizeof(std::complex<double>));
  h5_dcomplex.insertMember(REAL, 0, H5::PredType::NATIVE_DOUBLE);
  h5_dcomplex.insertMember(IMAG, sizeof(double), H5::PredType::NATIVE_DOUBLE);
  return h5_dcomplex;
}();

static void ReadAllElements(
    Eigen::Tensor<std::complex<double>, 4, Eigen::RowMajor>& coefficients,
    const H5::DataSet& dataset, const std::vector<unsigned int>& shape) {
  coefficients.resize(shape[0], shape[1], shape[2], shape[3]);
  dataset.read(coefficients.data(), kH5Dcomplex);
}

void ReadOneElement(
    Eigen::Tensor<std::complex<double>, 4, Eigen::RowMajor>& coefficients,
    const H5::DataSet& dataset, const std::vector<unsigned>& shape,
    unsigned index) {
  static constexpr size_t kRank = 4;

  // Define the part of the coefficients to read.
  const std::array<hsize_t, kRank> kOffset = {0, 0, index, 0};
  const std::array<hsize_t, kRank> kCount = {shape[0], shape[1], 1, shape[3]};
  static constexpr std::array<hsize_t, kRank> kStride = {1, 1, 1, 1};
  static constexpr std::array<hsize_t, kRank> kBlock = {1, 1, 1, 1};

  H5::DataSpace memspace{kRank, kCount.data()};
  H5::DataSpace dataspace = dataset.getSpace();
  dataspace.selectHyperslab(H5S_SELECT_SET, kCount.data(), kOffset.data(),
                            kStride.data(), kBlock.data());

  // TODO AST-807 The exact mapping between the data-layout of HD5 and
  // Eigen-Tensors needs to be investigated so the elements can be copied more
  // efficiently.  (Ideally they would be directly read in the proper shape.)
  std::vector<std::complex<double>> buffer(shape[0] * shape[1] * 1 * shape[3]);
  dataset.read(buffer.data(), kH5Dcomplex, memspace, dataspace);
  auto iterator = buffer.begin();
  coefficients.resize(shape[0], shape[1], 1, shape[3]);
  for (size_t i = 0; i < shape[0]; ++i) {
    for (size_t j = 0; j < shape[1]; ++j) {
      for (size_t k = 0; k < shape[3]; ++k) {
        coefficients(i, j, 0, static_cast<long>(k)) = *iterator++;
      }
    }
  }
#if 0
  // TODO AST-807 remove this validation code.
  // This validates the data has been read correctly when compared with the
  // original read function.
  Eigen::Tensor<std::complex<double>, 4, Eigen::RowMajor> expected;
  ReadAllElements(expected, dataset, shape);
  for (size_t i = 0; i < shape[0]; ++i) {
    for (size_t j = 0; j < shape[1]; ++j) {
      for (size_t k = 0; k < shape[3]; ++k) {
        if (coefficients(i, j, 0, k) != expected(i, j, index, k)) {
          asm("int3");
        }
      }
    }
  }
#endif
}

LOBESElementResponse::LOBESElementResponse(const std::string& name,
                                           const Options& options) {
  const std::optional<AartfaacStation> aartfaac_station =
      GetAartfaacStation(name, AartfaacElements::kInner);

  const std::filesystem::path search_path =
      options.coeff_path.empty() ? GetPath("lobes")
                                 : std::filesystem::path(options.coeff_path);
  const std::string_view station_name =
      aartfaac_station ? aartfaac_station->station : name;
  const std::string station_file = "LOBES_" + std::string(station_name) + ".h5";
  const std::filesystem::path coeff_file_path = search_path / station_file;
  H5::H5File h5file;

  if (!std::filesystem::exists(coeff_file_path)) {
    throw std::runtime_error("LOBES coeffcients file: " +
                             coeff_file_path.string() + " does not exists");
  }

  try {
    h5file.openFile(coeff_file_path.c_str(), H5F_ACC_RDONLY);
  } catch (const H5::FileIException& e) {
    throw std::runtime_error("Could not open LOBES coeffcients file: " +
                             coeff_file_path.string());
  }

  H5::DataSet dataset = h5file.openDataSet("coefficients");
  H5::DataSpace dataspace = dataset.getSpace();
  int nr_elements = dataspace.getSimpleExtentNpoints();

  // Get the number of dimensions in the dataspace.
  int ndims_coefficients = dataspace.getSimpleExtentNdims();

  // Get the dimension size of each dimension in the dataspace and display them.
  std::vector<hsize_t> dims_coefficients(ndims_coefficients);
  dataspace.getSimpleExtentDims(dims_coefficients.data(), nullptr);
  const std::vector<unsigned int> coefficients_shape(dims_coefficients.begin(),
                                                     dims_coefficients.end());

  if (aartfaac_station) {
    std::stringstream sstr;
    ReadOneElement(coefficients_, dataset, coefficients_shape,
                   aartfaac_station->element);
  } else {
    ReadAllElements(coefficients_, dataset, coefficients_shape);
  }

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

std::shared_ptr<ElementResponse> LOBESElementResponse::FixateDirection(
    const vector3r_t& direction) const {
  const vector2r_t thetaphi = cart2thetaphi(direction);

  return std::make_shared<LobesElementResponseFixedDirection>(
      std::static_pointer_cast<const LOBESElementResponse>(shared_from_this()),
      ComputeBaseFunctions(thetaphi[0], thetaphi[1]));
}

LOBESElementResponse::BaseFunctions LOBESElementResponse::ComputeBaseFunctions(
    double theta, double phi) const {
  BaseFunctions base_functions(nms_.size() * 2, 0.0);

  for (size_t i = 0; i < nms_.size(); ++i) {
    const nms_t& nms = nms_[i];
    std::complex<double> q2;
    std::complex<double> q3;
    std::tie(q2, q3) =
        everybeam::common::F4far_new(nms.s, nms.m, nms.n, theta, phi);
    base_functions[i * 2 + 0] = q2;
    base_functions[i * 2 + 1] = q3;
  }

  return base_functions;
}

aocommon::MC2x2 LOBESElementResponse::Response(int element_id, double frequency,
                                               double theta, double phi) const {
  // Clip directions below the horizon.
  if (theta >= M_PI_2) {
    return aocommon::MC2x2::Zero();
  }

  return Response(ComputeBaseFunctions(theta, phi), element_id, frequency);
}

aocommon::MC2x2 LOBESElementResponse::Response(
    const BaseFunctions& base_functions, int element_id,
    double frequency) const {
  const int frequency_index = FindFrequencyIndex(frequency);
  aocommon::MC2x2 response = aocommon::MC2x2::Zero();

  for (size_t i = 0; i < base_functions.size() / 2; ++i) {
    const std::complex<double> q2 = base_functions[i * 2 + 0];
    const std::complex<double> q3 = base_functions[i * 2 + 1];
    response[0] += q2 * coefficients_(0, frequency_index, element_id, i);  // xx
    response[1] += q3 * coefficients_(0, frequency_index, element_id, i);  // xy
    response[2] += q2 * coefficients_(1, frequency_index, element_id, i);  // yx
    response[3] += q3 * coefficients_(1, frequency_index, element_id, i);  // yy
  }

  return response;
}

std::shared_ptr<const LOBESElementResponse> LOBESElementResponse::GetInstance(
    const std::string& name, const Options& options) {
  // Using a single LOBESElementResponse object for each name reduces memory
  // usage since the coefficients are only loaded once.
  // Using weak pointers in this map ensures that LOBESElementResponse objects,
  // are deleted when they are no longer used, which saves memory.
  static std::map<std::string, std::weak_ptr<const LOBESElementResponse>>
      name_response_map;
  std::shared_ptr<const LOBESElementResponse> instance;

  auto entry = name_response_map.find(name);
  if (entry == name_response_map.end()) {
    instance = std::make_shared<const LOBESElementResponse>(name, options);
    name_response_map.insert({name, instance});
  } else {
    instance = entry->second.lock();
    if (!instance) {
      instance = std::make_shared<const LOBESElementResponse>(name, options);
      entry->second = instance;
    }
  }
  return instance;
}

}  // namespace everybeam
