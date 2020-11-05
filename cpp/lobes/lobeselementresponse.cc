#include <cmath>

#include "config.h"
#include "lobeselementresponse.h"

#include <map>
#include <tuple>
#include <complex>
#include <H5Cpp.h>
#include <iostream>

#include <boost/math/special_functions/legendre.hpp>

// Anonymous namespace for evaluating the base functions
namespace {
double P(int m, int n, double x) {
  double result = boost::math::legendre_p(n, std::abs(m), x);
  if (m < 0) {
    int phase = ((-m) % 2 == 0) ? 1 : -1;
    result *= phase * std::tgamma(n + m + 1) / std::tgamma(n - m + 1);
  }

  return result;
}

double Pacc(int m, int n, double x) {
  return (-(n + m) * (n - m + 1.0) * sqrt(1.0 - x * x) * P(m - 1, n, x) -
          m * x * P(m, n, x)) /
         (x * x - 1.0);
}

std::pair<std::complex<double>, std::complex<double>> F4far_new(int s, int m,
                                                                int n,
                                                                double theta,
                                                                double phi) {
  double C;
  if (m) {
    C = std::sqrt(60.0) * 1.0 / std::sqrt(n * (n + 1.0)) *
        std::pow(-m / std::abs(m), m);
  } else {
    C = std::sqrt(60.0) * 1.0 / std::sqrt(n * (n + 1.0));
  }

  std::complex<double> q2;
  std::complex<double> q3;

  // From cpp >= cpp14, complex literals can be used
  constexpr std::complex<double> i_neg = {0.0, -1.0};
  constexpr std::complex<double> i_pos = {0.0, 1.0};

  auto cos_theta = cos(theta);
  auto sin_theta = sin(theta);
  auto P_cos_theta = P(std::abs(m), n, cos_theta);
  auto Pacc_cos_theta = Pacc(std::abs(m), n, cos_theta);
  auto exp_i_m_phi = exp(i_pos * double(m) * phi);

  if (s == 1) {
    q2 = C * std::pow(i_neg, -n - 1) * i_pos * double(m) /
         (sin_theta)*std::sqrt((2. * n + 1) / 2.0 *
                               std::tgamma(n - std::abs(m) + 1) /
                               std::tgamma(n + std::abs(m) + 1)) *
         P_cos_theta * exp_i_m_phi;

    q3 = C * std::pow(i_neg, -n - 1) *
         std::sqrt((2. * n + 1) / 2.0 * std::tgamma(n - abs(m) + 1) /
                   std::tgamma(n + abs(m) + 1)) *
         Pacc_cos_theta * sin_theta * exp_i_m_phi;
  } else if (s == 2) {
    q2 = -C * std::pow(i_neg, -n) *
         std::sqrt((2. * n + 1) / 2.0 * std::tgamma(n - abs(m) + 1) /
                   std::tgamma(n + abs(m) + 1)) *
         Pacc_cos_theta * sin_theta * exp_i_m_phi;

    q3 = C * std::pow(i_neg, -n) * i_pos * double(m) / sin_theta *
         std::sqrt((2. * n + 1) / 2.0 * std::tgamma(n - abs(m) + 1) /
                   std::tgamma(n + abs(m) + 1)) *
         P_cos_theta * exp_i_m_phi;
  }

  return std::make_pair(q2, q3);
}
}  // namespace

namespace everybeam {
LOBESElementResponse::LOBESElementResponse(std::string name) {
  std::string file_name = std::string("LOBES_") + name + std::string(".h5");
  std::string path_name = GetPath(file_name.c_str());
  H5::H5File h5file;

  try {
    h5file.openFile(path_name.c_str(), H5F_ACC_RDONLY);
  } catch (const H5::FileIException& e) {
    throw std::runtime_error("Could not open LOBES coeffcients file: " +
                             path_name);
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
  dataspace.getSimpleExtentDims(dims_coefficients.data(), NULL);
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
  dataspace.getSimpleExtentDims(dims_nms.data(), NULL);

  nms_.resize(dims_nms[0]);
  dataset.read(nms_.data(), H5::PredType::NATIVE_INT);
}

LOBESElementResponse::BaseFunctions LOBESElementResponse::ComputeBaseFunctions(
    double theta, double phi) const {
  LOBESElementResponse::BaseFunctions base_functions(nms_.size(), 2);
  base_functions.setZero();

  // Avoid singularities due to theta = 0
  if (std::abs(theta) < 1e-6) {
    // TODO: should be handled by a proper warning
    std::cout << "LOBES is singular for zenith angle = 0rad. Will round "
              << theta << " to a value of 1e-6" << std::endl;
    theta = 1e-6;
  }

  for (size_t i = 0; i < nms_.size(); ++i) {
    auto nms = nms_[i];
    std::complex<double> q2, q3;
    std::tie(q2, q3) = F4far_new(nms.s, nms.m, nms.n, theta, phi);
    base_functions(i, 0) = q2;
    base_functions(i, 1) = q3;
  }
  return base_functions;
}

void LOBESElementResponse::Response(
    int element_id, double freq, double theta, double phi,
    std::complex<double> (&response)[2][2]) const {
  // Initialize the response to zero.
  response[0][0] = 0.0;
  response[0][1] = 0.0;
  response[1][0] = 0.0;
  response[1][1] = 0.0;

  // Clip directions below the horizon.
  if (theta >= M_PI_2) {
    return;
  }

  // Fill basefunctions if not yet set. Set clear_basefunctions to false
  // to disable the caching of the basefunctions
  bool clear_basefunctions = false;
  if (basefunctions_.rows() == 0) {
    clear_basefunctions = true;
    basefunctions_ = ComputeBaseFunctions(theta, phi);
  }

  int freq_idx = FindFrequencyIdx(freq);
  std::complex<double> xx = {0}, xy = {0}, yx = {0}, yy = {0};

  int nr_rows = basefunctions_.rows();
  if (nr_rows == 0) {
    throw std::runtime_error(
        "Number of rows in basefunctions_ member is 0. Did you run "
        "SetFieldQuantities?");
  }

  for (int i = 0; i < nr_rows; ++i) {
    std::complex<double> q2 = basefunctions_(i, 0);
    std::complex<double> q3 = basefunctions_(i, 1);
    xx += q2 * coefficients_(0, freq_idx, element_id, i);
    xy += q3 * coefficients_(0, freq_idx, element_id, i);
    yx += q2 * coefficients_(1, freq_idx, element_id, i);
    yy += q3 * coefficients_(1, freq_idx, element_id, i);
  }

  response[0][0] = xx;
  response[0][1] = xy;
  response[1][0] = yx;
  response[1][1] = yy;

  if (clear_basefunctions) {
    // Do a destructive resize
    basefunctions_.resize(0, 2);
  }
}

std::shared_ptr<LOBESElementResponse> LOBESElementResponse::GetInstance(
    std::string name) {
  static std::map<std::string, std::shared_ptr<LOBESElementResponse>>
      name_response_map;

  auto entry = name_response_map.find(name);
  if (entry == name_response_map.end()) {
    entry = name_response_map.insert(
        entry, {name, std::make_shared<LOBESElementResponse>(name)});
  }
  return entry->second;
}

std::string LOBESElementResponse::GetPath(const char* filename) const {
  std::stringstream ss;
  ss << EVERYBEAM_DATA_DIR << "/lobes/";
  ss << filename;
  return ss.str();
}

}  // namespace everybeam
