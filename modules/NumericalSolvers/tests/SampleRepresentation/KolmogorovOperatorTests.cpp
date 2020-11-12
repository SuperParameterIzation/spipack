#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/NumericalSolvers/SampleRepresentation/KolmogorovOperator.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

class KolmogorovOperatorTests : public::testing::Test {
protected:

  /// Set up information to test the sample representation
  virtual void SetUp() override {
    // create a standard Gaussian random variable
    rv = std::make_shared<Gaussian>(dim)->AsVariable();

    // options for the nearest neighbor search
    YAML::Node nnOptions;
    nnOptions["NumSamples"] = n;
    nnOptions["Stride"] = n/5;
    nnOptions["NumThreads"] = omp_get_max_threads();

    // set the kernel options
    YAML::Node kernelOptions;
    kernelOptions["Kernel"] = "ExponentialKernel";

    // set the options for the Kolmogorov operator
    options["NearestNeighbors"] = nnOptions;
    options["NumNearestNeighbors"] = nneighs;
    options["KernelOptions"] = kernelOptions;
  }

  /// Make sure everything is what we expect
  virtual void TearDown() override {
    EXPECT_EQ(kolOperator->NumNearestNeighbors(), nneighs);
    EXPECT_EQ(kolOperator->NumSamples(), n);
  }

  /// Create the sample representation from samples
  /**
    \return The sample collection used to create the Kolmogorov operator
  */
  inline std::shared_ptr<SampleCollection> CreateFromSamples() {
    // add random samples into a sample collection
    auto samples = std::make_shared<SampleCollection>();
    for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

    // create the Kolmogorov operator
    kolOperator = std::make_shared<KolmogorovOperator>(samples, options);

    // return the samples
    return samples;
  }

  /// The dimension of state spaces
  const unsigned int dim = 4;

  /// The number of samples
  std::size_t n = 1000;

  /// The number of nearest neighbors
  const std::size_t nneighs = 15;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The options for the graph Laplacian
  YAML::Node options;

  /// The Kolmogorov operator---use a pointer here so we can initalize it as null
  std::shared_ptr<KolmogorovOperator> kolOperator;
};

TEST_F(KolmogorovOperatorTests, RandomVariableConstruction) {
  // create the Kolmogorov operator
  kolOperator = std::make_shared<KolmogorovOperator>(rv, options);
}

TEST_F(KolmogorovOperatorTests, SampleCollectionConstruction) {
  // create the Kolmogorov operator from samples
  auto samples = CreateFromSamples();

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-kolOperator->Point(i)).norm(), 0.0, 1.0e-10);
  }
}
