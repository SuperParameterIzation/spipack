#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/NumericalSolvers/GraphLaplacian/SampleRepresentation.hpp"

using namespace muq::Modeling;
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

  /// The dimension of state spaces
  const unsigned int dim = 4;

  /// The number of samples
  const std::size_t n = 2000;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The options for the graph Laplacian
  YAML::Node options;

private:
};

TEST_F(SampleRepresentationTests, RandomVariableConstruction) {

}
