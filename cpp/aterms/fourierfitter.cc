#include "fourierfitter.h"

#include <xtensor/xadapt.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor-fftw/basic.hpp>
#include <xtensor/xvectorize.hpp>

namespace everybeam {
namespace aterms {

FourierFitter::FourierFitter(
    std::size_t subgrid_size, std::size_t support,
    const std::vector<std::pair<float, float>>& directions)
    : subgrid_size_(subgrid_size), support_(support) {
  // The solution to the least squares fitting problem solved in Evaluate() is
  // a multiplication by a matrix inverse.
  // This matrix depends on the directions, which are known here.
  // The inverse matrix can be pre-computed here once and applied to the data
  // in all subsequent calls to Evaluate()
  const size_t nr_unknowns = support * support;
  const size_t nr_directions = directions.size();

  xt::xtensor<std::complex<float>, 2>::shape_type shape = {nr_directions,
                                                           nr_unknowns};

  xt::xtensor<std::complex<float>, 2> a(shape);

  for (size_t j = 0; j < nr_directions; ++j) {
    auto [x, y] = directions[j];
    for (size_t k = 0; k < nr_unknowns; ++k) {
      const int ii = int(k) / support - support / 2;
      const int jj = int(k) % support - support / 2;
      a[{j, k}] = std::exp(std::complex<float>(
          0, -2 * M_PI * (ii * x + jj * y) / float(subgrid_size)));
    }
  }

  a_inv_ = xt::linalg::pinv(a);
}

void FourierFitter::Evaluate(const std::vector<std::complex<float>>& solutions,
                             std::complex<float>* buffer) const {
  // Solve matrix equation Ax = y
  // where y are the calibration solutions.
  // Matrix A is the same for all calls to Evaluate() and
  // therefore the inverse of matrix A is precomputed and stored in member
  // a_inv_

  // Create a vector (N x 1 tensor) from a std::vector
  // xt::adapt converts from std containers to xtensors
  // The new xt::view is of all pre-existing dimensions plus a new additional
  // length-1 axis
  xt::xtensor<std::complex<float>, 2> y =
      xt::view(xt::adapt(solutions), xt::all(), xt::newaxis());

  // Compute matrix vector product A^{-1} y
  xt::xtensor<std::complex<float>, 2> x = xt::linalg::dot(a_inv_, y);

  // Create an empty screen in the inverse Fourier transformed domain
  xt::xarray<std::complex<float>> screen_ift =
      xt::zeros<std::complex<float>>({subgrid_size_, subgrid_size_});

  // Fill the centre of the inverse Fourier transform of the screen
  const size_t nr_unknowns = support_ * support_;
  for (size_t k = 0; k < nr_unknowns; ++k) {
    const int ii = int(k) / support_ - support_ / 2;
    const int jj = int(k) % support_ - support_ / 2;
    screen_ift.periodic(ii, jj) = x(k, 0);
  }

  // Transform the screen to image space
  xt::xarray<std::complex<float>> screen = xt::fftw::fft2(screen_ift);

  // Write the result into the buffer
  const std::size_t size = subgrid_size_ * subgrid_size_;
  std::vector<std::size_t> shape = {subgrid_size_, subgrid_size_};
  xt::adapt(buffer, size, xt::no_ownership(), shape) = screen;
}

}  // namespace aterms
}  // namespace everybeam
