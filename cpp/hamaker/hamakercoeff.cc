// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "hamakercoeff.h"

namespace everybeam {
H5::CompType GetComplexDoubleType() {
  H5::CompType complex_type(sizeof(std::complex<double>));
  complex_type.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
  complex_type.insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);
  return complex_type;
}

size_t HamakerCoefficients::GetIndex(unsigned int n, unsigned int t,
                                     unsigned int f) const {
  return (n * nPowerTheta_ + t) * nPowerFreq_ * nInner_ + f * nInner_;
}

// Constructor for reading coeff from file
HamakerCoefficients::HamakerCoefficients(const std::string& filename) {
  ReadCoefficients(filename);
}

// Constructor for writing coeff to file
HamakerCoefficients::HamakerCoefficients(const double freq_center,
                                         const double freq_range,
                                         const unsigned int nHarmonics,
                                         const unsigned int nPowerTheta,
                                         const unsigned int nPowerFreq)
    : freq_center_(freq_center),
      freq_range_(freq_range),
      nHarmonics_(nHarmonics),
      nPowerTheta_(nPowerTheta),
      nPowerFreq_(nPowerFreq),
      coeff_(GetNrCoefficients()) {}

size_t HamakerCoefficients::GetNrCoefficients() const {
  return nHarmonics_ * nPowerTheta_ * nPowerFreq_ * nInner_;
}

void HamakerCoefficients::SetCoefficients(
    unsigned int n, unsigned int t, unsigned int f,
    std::pair<std::complex<double>, std::complex<double>> value) {
  size_t index = GetIndex(n, t, f);
  coeff_[index] = value.first;
  coeff_[index + 1] = value.second;
}

void HamakerCoefficients::SetCoefficients(const std::complex<double>* coeff) {
  memcpy(coeff_.data(), coeff, coeff_.size() * sizeof(std::complex<double>));
}

void HamakerCoefficients::SetCoefficients(
    const std::vector<std::complex<double>> coeff) {
  assert(coeff.size() == coeff_.size());
  std::copy(coeff.begin(), coeff.end(), coeff_.begin());
}

void HamakerCoefficients::GetCoefficient(
    unsigned int n, unsigned int t, unsigned int f,
    std::pair<std::complex<double>, std::complex<double>>& value) const {
  size_t index = GetIndex(n, t, f);
  value.first = coeff_[index];
  value.second = coeff_[index + 1];
}

void HamakerCoefficients::ReadCoefficients(const std::string& filename) {
  // Open file
  H5::H5File file(filename, H5F_ACC_RDONLY);

  // Read dataset
  H5::DataSet dataset = file.openDataSet(dataset_name_);

  // Open attribute and read its contents
  H5::Attribute freq_center_attr = dataset.openAttribute("freq_center");
  H5::Attribute freq_range_attr = dataset.openAttribute("freq_range");
  freq_center_attr.read(H5::PredType::NATIVE_DOUBLE, &freq_center_);
  freq_range_attr.read(H5::PredType::NATIVE_DOUBLE, &freq_range_);

  // Read dataset dimensions
  H5::DataSpace dataspace = dataset.getSpace();
  unsigned int rank = dataspace.getSimpleExtentNdims();
  assert(rank == dataset_rank_);
  std::vector<hsize_t> dims(rank);
  dataspace.getSimpleExtentDims(dims.data(), nullptr);
  nHarmonics_ = dims[0];
  nPowerTheta_ = dims[1];
  nPowerFreq_ = dims[2];

  // Read coeffs
  coeff_.resize(GetNrCoefficients());
  H5::DataType data_type = dataset.getDataType();
  assert(data_type.getSize() == sizeof(std::complex<double>));
  dataset.read(coeff_.data(), data_type, dataspace);
}

void HamakerCoefficients::WriteCoefficients(const std::string& filename) {
  // Open file
  H5::H5File file(filename, H5F_ACC_TRUNC);

  // Create dataspace
  const unsigned int rank = 4;
  hsize_t dims[rank] = {nHarmonics_, nPowerTheta_, nPowerFreq_, nInner_};
  H5::DataSpace coeff_dataspace(rank, dims);

  // Create complex type
  H5::CompType complex_type = GetComplexDoubleType();

  // Write dataset
  H5::DataSet dataset =
      file.createDataSet("coeff", complex_type, coeff_dataspace);
  dataset.write(coeff_.data(), complex_type);

  // Create new dataspace for attribute
  H5::DataSpace attr_dataspace(H5S_SCALAR);

  // Write frequency center attribute
  H5::Attribute freq_center_attr = dataset.createAttribute(
      "freq_center", H5::PredType::NATIVE_DOUBLE, attr_dataspace);
  freq_center_attr.write(H5::PredType::NATIVE_DOUBLE, &freq_center_);

  // Write frequency range attribute
  H5::Attribute freq_range_attr = dataset.createAttribute(
      "freq_range", H5::PredType::NATIVE_DOUBLE, attr_dataspace);
  freq_range_attr.write(H5::PredType::NATIVE_DOUBLE, &freq_range_);

  file.flush(H5F_SCOPE_LOCAL);
}

void HamakerCoefficients::PrintCoefficients() {
  std::pair<std::complex<double>, std::complex<double>> coeff;
  for (unsigned int h = 0; h < nHarmonics_; h++) {
    for (unsigned int t = 0; t < nPowerTheta_; t++) {
      for (unsigned int f = 0; f < nPowerFreq_; f++) {
        GetCoefficient(h, t, f, coeff);
        std::cout << coeff.first << ", " << coeff.second << std::endl;
      }
    }
  }
  std::cout << std::endl;
}
}  // namespace everybeam
