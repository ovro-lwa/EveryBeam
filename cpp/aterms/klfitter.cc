#include "klfitter.h"

#include <xtensor/xadapt.hpp>
#include <xtensor/xio.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor-fftw/basic.hpp>
#include <xtensor/xvectorize.hpp>

namespace everybeam {
namespace aterms {

KlFitter::KlFitter(std::size_t subgrid_size, int order,
                   const std::vector<std::pair<float, float>>& directions)
    : subgrid_size_(subgrid_size) {
  std::vector<size_t> shape = {directions.size(), 2};
  size_t size = directions.size() * 2;

  auto directions_tensor =
      xt::adapt(reinterpret_cast<const float*>(directions.data()), size,
                xt::no_ownership(), shape);

  // Compute the difference vector per pair of directions
  // Use two different views and broadcasting to expand from
  // all directions to all pairs of directions.
  auto directions_view1 =
      xt::view(directions_tensor, xt::all(), xt::newaxis(), xt::all());
  auto directions_view2 =
      xt::view(directions_tensor, xt::newaxis(), xt::all(), xt::all());
  auto difference = directions_view1 - directions_view2;

  auto distance_squared = xt::sum(difference * difference, {2});

  // Exponent of the power-law distriution for Kolmogorov turbulence
  const float beta = 5. / 3;
  // Use beta divided by two here, because the distance is already squared
  auto corr_matrix = xt::pow(distance_squared, beta / 2);
  auto [u, s, v] = xt::linalg::svd(corr_matrix);
  auto basis_vectors = xt::view(u, xt::all(), xt::range(0, order));
  auto projection_matrix =
      xt::linalg::dot(basis_vectors, xt::transpose(basis_vectors));
  auto corr_matrix_inv = xt::linalg::pinv(corr_matrix, 1e-3);

  // Positions of all pixels in the screen
  auto pixel_positions = xt::stack(
      xt::meshgrid(xt::arange(subgrid_size), xt::arange(subgrid_size)), 2);

  // Difference vectors between all pixel-direction(source) pairs
  auto difference_pixel_source = xt::view(pixel_positions, xt::all(), xt::all(),
                                          xt::newaxis(), xt::all()) -
                                 xt::view(directions_tensor, xt::newaxis(),
                                          xt::newaxis(), xt::all(), xt::all());

  auto distance_squared_pixel_source =
      xt::sum(difference_pixel_source * difference_pixel_source, -1);

  // Cross correlations between pixels and sources
  auto cross_corr_matrix_pixel_source =
      xt::pow(distance_squared_pixel_source, beta / 2);

  fitting_matrix_ = xt::linalg::tensordot(
      cross_corr_matrix_pixel_source,
      xt::linalg::dot(corr_matrix_inv, projection_matrix), 1);
}

void KlFitter::Evaluate(const std::vector<float>& solutions,
                        float* buffer) const {
  // Convert 'buffer' and 'solutions' into xtensor objects and compute the
  // matrix-vector product.
  const std::size_t size = subgrid_size_ * subgrid_size_;
  const std::array<std::size_t, 2> shape = {subgrid_size_, subgrid_size_};
  xt::adapt(buffer, size, xt::no_ownership(), shape) =
      xt::linalg::tensordot(fitting_matrix_, xt::adapt(solutions), 1);
}

}  // namespace aterms
}  // namespace everybeam
