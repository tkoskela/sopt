#include <algorithm>
#include <exception>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <ctime>

#include <sopt/imaging_padmm.h>
#include <sopt/logging.h>
#include <sopt/maths.h>
#include <sopt/positive_quadrant.h>
#include <sopt/relative_variation.h>
#include <sopt/reweighted.h>
#include <sopt/sampling.h>
#include <sopt/types.h>
#include <sopt/utilities.h>
#include <sopt/wavelets.h>
#include <sopt/wavelets/sara.h>
// This header is not part of the installed sopt interface
// It is only present in tests
#include <tools_for_tests/directories.h>
#include <tools_for_tests/tiffwrappers.h>

// \min_{x} ||\Psi^Tx||_1 \quad \mbox{s.t.} \quad ||y - Ax||_2 < \epsilon and x \geq 0
int main(int argc, char const **argv) {
  // Some typedefs for simplicity
  typedef double Scalar;
  // Column vector - linear algebra - A * x is a matrix-vector multiplication
  // type expected by ProximalADMM
  typedef sopt::Vector<Scalar> Vector;
  // Matrix - linear algebra - A * x is a matrix-vector multiplication
  // type expected by ProximalADMM
  typedef sopt::Matrix<Scalar> Matrix;
  // Image - 2D array - A * x is a coefficient-wise multiplication
  // Type expected by wavelets and image write/read functions
  typedef sopt::Image<Scalar> Image;

  std::string const input = argc >= 2 ? argv[1] : "cameraman256";
  std::string const output = argc == 3 ? argv[2] : "none";
  if (argc > 3) {
    std::cout << "Usage:\n"
                 "$ "
              << argv[0]
              << " [input [output]]\n\n"
                 "- input: path to the image to clean (or name of standard SOPT image)\n"
                 "- output: filename pattern for output image\n";
    exit(0);
  }
  // Set up random numbers for C and C++
  auto const seed = std::time(0);
  std::srand((unsigned int)seed);
  std::mt19937 mersenne(std::time(0));

  // Initializes and sets logger (if compiled with logging)
  // See set_level function for levels.
  sopt::logging::initialize();

  SOPT_MEDIUM_LOG("Read input file {}", input);
  Image const image = sopt::notinstalled::read_standard_tiff(input);

  SOPT_MEDIUM_LOG("Initializing sensing operator");
  sopt::t_uint nmeasure = 0.33 * image.size();
  auto const sampling =
      sopt::linear_transform<Scalar>(sopt::Sampling(image.size(), nmeasure, mersenne));

  SOPT_MEDIUM_LOG("Initializing wavelets");
  sopt::wavelets::SARA const sara{std::make_tuple(std::string{"DB3"}, 1u),
                                  std::make_tuple(std::string{"DB1"}, 2u),
                                  std::make_tuple(std::string{"DB1"}, 3u)};
  auto const psi = sopt::linear_transform<Scalar>(sara, image.rows(), image.cols());

  SOPT_MEDIUM_LOG("Computing proximal-ADMM parameters");
  Vector const y0 = sampling * Vector::Map(image.data(), image.size());
  auto const snr = 30.0;
  auto const sigma = y0.stableNorm() / std::sqrt(y0.size()) * std::pow(10.0, -(snr / 20.0));
  auto const epsilon = std::sqrt(nmeasure + 2 * std::sqrt(y0.size())) * sigma;

  SOPT_MEDIUM_LOG("Create dirty vector");
  std::normal_distribution<> gaussian_dist(0, sigma);
  Vector y(y0.size());
  for (sopt::t_int i = 0; i < y0.size(); i++) y(i) = y0(i) + gaussian_dist(mersenne);
  // Write dirty imagte to file
  if (output != "none") {
    Vector const dirty = sampling.adjoint() * y;
    sopt::utilities::write_tiff(Matrix::Map(dirty.data(), image.rows(), image.cols()),
                                "dirty_" + output + ".tiff");
  }

  SOPT_MEDIUM_LOG("Creating proximal-ADMM Functor");
  auto const padmm = sopt::algorithm::ImagingProximalADMM<Scalar>(y)
                         .itermax(500)
                         .gamma(1e-1)
                         .relative_variation(5e-4)
                         .l2ball_proximal_epsilon(epsilon)
                         .tight_frame(false)
                         .l1_proximal_tolerance(1e-2)
                         .l1_proximal_nu(1)
                         .l1_proximal_itermax(50)
                         .l1_proximal_positivity_constraint(true)
                         .l1_proximal_real_constraint(true)
                         .residual_convergence(epsilon * 1.001)
                         .lagrange_update_scale(0.9)
                         .nu(1e0)
                         .Psi(psi)
                         .Phi(sampling);

  SOPT_MEDIUM_LOG("Creating the reweighted algorithm");
  // positive_quadrant projects the result of PADMM on the positive quadrant.
  // This follows the reweighted algorithm for SDMM
  auto const min_delta = sigma * std::sqrt(y.size()) / std::sqrt(8 * image.size());
  auto const reweighted =
      sopt::algorithm::reweighted(padmm).itermax(5).min_delta(min_delta).is_converged(
          sopt::RelativeVariation<Scalar>(1e-3));

  SOPT_MEDIUM_LOG("Starting proximal-ADMM");
  // Alternatively, padmm can be called with a tuple (x, residual) as argument
  // Here, we default to (Φ^Ty/ν, ΦΦ^Ty/ν - y)
  auto const diagnostic = reweighted();
  SOPT_MEDIUM_LOG("proximal-ADMM returned {}", diagnostic.good);

  SOPT_MEDIUM_LOG("SOPT-proximal-ADMM converged in {} iterations", diagnostic.niters);
  if (output != "none")
    sopt::utilities::write_tiff(Matrix::Map(diagnostic.algo.x.data(), image.rows(), image.cols()),
                                output + ".tiff");

  return 0;
}
