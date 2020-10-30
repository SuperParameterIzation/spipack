#include <gtest/gtest.h>

#include "spipack/NumericalSolvers/GraphLaplacian.hpp"

#include <MUQ/Modeling/Distributions/Gaussian.h>

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

class GraphLaplacianTests : public::testing::Test {
public:
  /// Set up information to test the graph Laplacian
  void SetUp() {
    // create a standard Gaussian random variable
    rv = std::make_shared<Gaussian>(dim)->AsVariable();

    // set the options for the graph laplacian
    options["MaxLeaf"] = maxLeaf;
  }

  /*/// Make sure everything is constructed correctly
  virtual void TearDown() override {
    // make sure the graph laplacian has enough samples
    EXPECT_EQ(laplacian->NumSamples(), n);

    // make the the kd tree max leaf is the value we set
    EXPECT_EQ(laplacian->KDTreeMaxLeaf(), maxLeaf);
  }*/

protected:
  /// The dimension of state space
  inline static const unsigned int dim = 4;

  /// The number of samples
  const size_t n = 1000;

  /// The max leaf size for the kd tree
  const size_t maxLeaf = 15;

  /// The options for the graph Laplacian
  YAML::Node options;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The graph Laplacian---use a pointer here so we can initalize it as null
  std::shared_ptr<GraphLaplacian> laplacian;
};

//TEST_F(GraphLaplacianTests, RandomVariableConstruction) {
  // create the options for the graph laplacian
  //options["NumSamples"] = n;

  // create the graph laplacian
  //laplacian = std::make_shared<GraphLaplacian>(rv, options);
//}

TEST(TESTGraphLaplacianTests, SampleCollectionConstruction) {
  /// The dimension of state space
  static const unsigned int dim = 4;

  /// The number of samples
  const size_t n = 1000;

  /// The max leaf size for the kd tree
  const size_t maxLeaf = 15;

  /// The options for the graph Laplacian
  YAML::Node options;

  // create a standard Gaussian random variable
  auto rv = std::make_shared<Gaussian>(dim)->AsVariable();

  // set the options for the graph laplacian
  options["MaxLeaf"] = maxLeaf;

  // add random samples into a sample collection
  auto samples = std::make_shared<SampleCollection>();
  for( size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

  // create the graph laplacian
  auto laplacian = std::make_shared<GraphLaplacian>(samples, options);

  // check to make sure the samples match
  /*for( size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-laplacian->Point(i)).norm(), 0.0, 1.0e-10);
  }*/
}
