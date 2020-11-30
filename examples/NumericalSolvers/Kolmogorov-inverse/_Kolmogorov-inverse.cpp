#include <iostream>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <spipack/NumericalSolvers/SampleRepresentation/KolmogorovOperator.hpp>

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace spi::NumericalSolvers;

int main(int argc, char **argv) {
  // the dimension of the problem
  const std::size_t dim = 2;

  // the output filename
  const std::string filename = "outputData.h5";

  // the number of samples
  const std::size_t n = 1000;

  // numerical parameters
  const std::size_t numNeighbors = 25;

  // the number of eigenvalues
  const std::size_t neigs = 100;
  const double eigensolverTol = 1.0e-8;

  // create a standard Gaussian random variable
  auto rv = std::make_shared<Gaussian>(dim)->AsVariable();

  // options for the nearest neighbor search
  YAML::Node nnOptions;
  nnOptions["NumSamples"] = n;
  nnOptions["Stride"] = n/5;
  nnOptions["NumThreads"] = omp_get_max_threads();

  // set the kernel options
  YAML::Node kernelOptions;
  kernelOptions["Kernel"] = "ExponentialKernel";

  // set the options for the density estimation
  YAML::Node densityOptions;
  densityOptions["NearestNeighbors"] = nnOptions;
  densityOptions["KernelOptions"] = kernelOptions;
  densityOptions["NumNearestNeighbors"] = numNeighbors;
  densityOptions["TruncationTolerance"] = -std::log(1.0e-2);
  densityOptions["ManifoldDimension"] = (double)dim;

  // set the options for the Kolmogorov operator
  YAML::Node options;
  options["NearestNeighbors"] = nnOptions;
  options["KernelOptions"] = kernelOptions;
  options["DensityOptions"] = densityOptions;
  options["OperatorParameter"] = 1.0;
  options["TruncationTolerance"] = -std::log(1.0e-2);
  options["ManifoldDimension"] = (double)dim;
  options["BandwidthExponent"] = -0.5;
  options["BandwidthParameter"] = 1.0e-1;
  options["NumEigenvalues"] = neigs;
  options["EigensolverTolerance"] = eigensolverTol;

  // create the Kolmogorov operator
  auto kolOperator = std::make_shared<KolmogorovOperator>(rv, options);

  // construct the kd-trees
  std::cout << "building kd trees ... " << std::flush;
  kolOperator->BuildKDTrees();
  std::cout << "done." << std::endl;

  // tune the bandwidth parameter
  std::cout << "tune the bandwidth parameter ... " << std::flush;
  kolOperator->TuneBandwidthParameter(true);
  std::cout << "done." << std::endl;

  // compute the eigendecomposition of Lhat
  std::cout << "computing eigendecomposition of Lhat ... " << std::flush;
  Eigen::VectorXd Sinv(n), lambda(neigs);
  Eigen::MatrixXd Uhat(n, neigs);
  kolOperator->ComputeEigendecomposition(Sinv, lambda, Uhat);
  std::cout << "done." << std::endl;
}
