#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/NumericalSolvers/GraphLaplacian/SampleRepresentation.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

class SampleRepresentationTests : public::testing::Test {
public:

  /// Set up information to test the graph Laplacian
  virtual void SetUp() override {
    // create a standard Gaussian random variable
    rv = std::make_shared<Gaussian>(dim)->AsVariable();

    // options for the nearest neighbor search
    YAML::Node nnOptions;
    nnOptions["NumSamples"] = n;
    nnOptions["Stride"] = n/5;

    // set the options for the graph laplacian
    options["NearestNeighbors"] = nnOptions;
  }

protected:

  /// Create the graph Laplacian from samples
  /**
    \return The sample collection used to create the graph Laplacian
  */
  inline std::shared_ptr<SampleCollection> CreateFromSamples() {
    // add random samples into a sample collection
    auto samples = std::make_shared<SampleCollection>();
    for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

    // create the graph laplacian
    representation = std::make_shared<SampleRepresentation>(samples, options);

    // return the samples
    return samples;
  }

  /// The dimension of state spaces
  const unsigned int dim = 4;

  /// The number of samples
  const std::size_t n = 2000;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The options for the graph Laplacian
  YAML::Node options;

  /// The sample representation---use a pointer here so we can initalize it as null
  std::shared_ptr<SampleRepresentation> representation;

private:
};

TEST_F(SampleRepresentationTests, RandomVariableConstruction) {
  // create the graph laplacian
  representation = std::make_shared<SampleRepresentation>(rv, options);
}

TEST_F(SampleRepresentationTests, SampleCollectionConstruction) {
  // create the graph laplacian from samples
  auto samples = CreateFromSamples();

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-representation->Point(i)).norm(), 0.0, 1.0e-10);
  }
}
