#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/Density.h>
#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

class DensityEstimationTests : public::testing::Test {
protected:

  /// Set up information to test the sample representation
  virtual void SetUp() override {
    // create a standard Gaussian random variable
    auto gauss = std::make_shared<Gaussian>(dim);
    rv = gauss->AsVariable();
    dens = gauss->AsDensity();

    // options for the nearest neighbor search
    YAML::Node nnOptions;
    nnOptions["NumSamples"] = n;
    nnOptions["Stride"] = n/5;
    nnOptions["NumThreads"] = omp_get_max_threads();

    // set the kernel options
    YAML::Node kernelOptions;
    kernelOptions["Kernel"] = "HatKernel";

    // set the options for the graph laplacian
    options["NearestNeighbors"] = nnOptions;
    options["NumNearestNeighbors"] = nneighs;
    options["KernelOptions"] = kernelOptions;
    options["BandwidthParameter"] = eps;
    options["ManifoldDimension"] = (double)dim;
  }

  /// Make sure everything is what we expect
  virtual void TearDown() override {
    EXPECT_EQ(density->NumNearestNeighbors(), nneighs);
    EXPECT_EQ(density->NumSamples(), n);
    EXPECT_DOUBLE_EQ(density->BandwidthParameter(), eps);
  }

  /// Create the density estimation from samples
  /**
    \return The sample collection used to create the graph Laplacian
  */
  inline std::shared_ptr<SampleCollection> CreateFromSamples() {
    // add random samples into a sample collection
    auto samples = std::make_shared<SampleCollection>();
    for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

    // create the graph laplacian
    density = std::make_shared<DensityEstimation>(samples, options);

    // return the samples
    return samples;
  }

  /// The dimension of state spaces
  const unsigned int dim = 1;

  /// The number of samples
  const std::size_t n = 10000;

  /// The number of nearest neighbors
  const std::size_t nneighs = 500;

  /// The bandwidth parameter \f$\epsilon\f$
  const double eps = 25.0;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The density associated with the random variable
  std::shared_ptr<Density> dens;

  /// The options for the graph Laplacian
  YAML::Node options;

  /// The sample representation---use a pointer here so we can initalize it as null
  std::shared_ptr<DensityEstimation> density;
};

TEST_F(DensityEstimationTests, RandomVariableConstruction) {
  // create the graph laplacian
  density = std::make_shared<DensityEstimation>(rv, options);
}

TEST_F(DensityEstimationTests, SampleCollectionConstruction) {
  // create the graph laplacian from samples
  auto samples = CreateFromSamples();

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-density->Point(i)).norm(), 0.0, 1.0e-10);
  }
}

/*
TEST_F(DensityEstimationTests, EstimateGaussian) {
  // create the graph laplacian from samples
  auto samples = CreateFromSamples();
  EXPECT_EQ(samples->size(), n);

  // construct the kd-trees
  density->BuildKDTrees();

  // create the tuning data
  DensityEstimation::TuningData tune;
  tune.bandwidthExponent = Eigen::VectorXd::LinSpaced(25, -3.0, 3.0);

  // estimate the density at each sample
  const Eigen::VectorXd densityEstimate = density->Estimate(tune);

  // the coefficient of the max. density point
  int coeff;
  double map = densityEstimate.maxCoeff(&coeff);

  //std::cout << coeff << std::endl;
  //std::cout << samples->at(coeff)->state[0].transpose() << std::endl;
  //std::cout << map << " " << std::exp(dens->LogDensity(samples->at(coeff)->state[0])) << std::endl;

  // check the estimate
  EXPECT_EQ(densityEstimate.size(), n);
}
*/
