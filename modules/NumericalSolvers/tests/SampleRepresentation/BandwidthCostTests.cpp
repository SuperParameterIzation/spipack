#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/NumericalSolvers/SampleRepresentation/BandwidthCost.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

TEST(BandwidthCostTests, ConstructCost) {
  // the parameter epsilon
  const double eps = std::log(10.0);

  // The dimension of state spaces
  const unsigned int dim = 2;

  // The number of samples
  std::size_t n = 1000;

  // The number of nearest neighbors
  const std::size_t nneighs = 15;

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

  // set the options for the sample representation
  YAML::Node options;
  options["NearestNeighbors"] = nnOptions;
  options["NumNearestNeighbors"] = nneighs;
  options["KernelOptions"] = kernelOptions;

  // add random samples into a sample collection
  auto samples = std::make_shared<SampleCollection>();
  for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

  // create a nearest neighbors object and compute the bandwidth
  NearestNeighbors nn(samples, options["NearestNeighbors"]);
  std::vector<std::vector<std::pair<std::size_t, double> > > neighbors;
  nn.BuildKDTrees();
  const Eigen::VectorXd bandwidth = (nn.SquaredBandwidth(nneighs, neighbors)).array().sqrt();

  // create the sample representation
  auto representation = std::make_shared<SampleRepresentation>(samples, options);
  representation->BuildKDTrees();

  // create the cost function
  BandwidthCost bandCost(bandwidth, representation);

  // compute the cost gradient using exact and FD
  const Eigen::VectorXd costGrad = bandCost.Gradient(0, Eigen::VectorXd::Constant(1, eps).eval(), Eigen::VectorXd::Ones(1).eval());
  const double gradFD = (bandCost.Cost(Eigen::VectorXd::Constant(1, eps+1.0e-8).eval())-bandCost.Cost(Eigen::VectorXd::Constant(1, eps).eval()))/1.0e-8;
  EXPECT_EQ(costGrad.size(), 1);
  EXPECT_NEAR(costGrad(0), gradFD, 1.0e-2);
}
