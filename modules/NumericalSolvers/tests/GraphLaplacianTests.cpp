#include <gtest/gtest.h>

#include "spipack/NumericalSolvers/GraphLaplacian.hpp"

#include <MUQ/Modeling/Distributions/Gaussian.h>

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

class GraphLaplacianTests : public::testing::Test {
public:
  /// Set up information to test the graph Laplacian
  virtual void SetUp() override {
    // create a standard Gaussian random variable
    rv = std::make_shared<Gaussian>(dim)->AsVariable();

    // set the options for the graph laplacian
    options["MaxLeaf"] = maxLeaf;
    options["NumSamples"] = n;
  }

  /// Make sure everything is constructed correctly
  virtual void TearDown() override {
    // make sure the graph laplacian has enough samples
    EXPECT_EQ(laplacian->NumSamples(), n);

    // make the the kd tree max leaf is the value we set
    EXPECT_EQ(laplacian->KDTreeMaxLeaf(), maxLeaf);
  }

protected:
  /// Create the graph Laplacian from samples
  /**
    \return The sample collection used to create the graph Laplacian
  */
  std::shared_ptr<SampleCollection> CreateFromSamples() {
    // add random samples into a sample collection
    auto samples = std::make_shared<SampleCollection>();
    for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

    // create the graph laplacian
    laplacian = std::make_shared<GraphLaplacian>(samples, options);

    // return the samples
    return samples;
  }

  /// The dimension of state spaces
  inline static const unsigned int dim = 4;

  /// The number of samples
  const std::size_t n = 10000;

  /// The max leaf size for the kd tree
  const std::size_t maxLeaf = 15;

  /// The options for the graph Laplacian
  YAML::Node options;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The graph Laplacian---use a pointer here so we can initalize it as null
  std::shared_ptr<GraphLaplacian> laplacian;
};

TEST_F(GraphLaplacianTests, RandomVariableConstruction) {
  // create the graph laplacian
  laplacian = std::make_shared<GraphLaplacian>(rv, options);
}

TEST_F(GraphLaplacianTests, SampleCollectionConstruction) {
  // create the graph laplacian from samples
  auto samples = CreateFromSamples();

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-laplacian->Point(i)).norm(), 0.0, 1.0e-10);
  }
}

TEST_F(GraphLaplacianTests, FindNearestNeighbors_Radius) {
  // create the graph laplacian from samples
  auto samples = CreateFromSamples();

  // build the kd-tree based on the samples
  laplacian->BuildKDTree();

  // choose a new random point from the distribution
  const Eigen::VectorXd x = rv->Sample();
  const double r = 0.5; // the radius

  // find all of the points within a choosen radius
  std::vector<std::pair<std::size_t, double> > neighbors;
  laplacian->FindNeighbors(x, r, neighbors);

  // make sure the index of the neighbors is valid and they are within the correct radius
  for( const auto& it : neighbors ) {
    EXPECT_TRUE(it.first<n);
    EXPECT_TRUE(it.second<r*r+1.0e-10);

    // make sure that we have found the correct neighbor
    const Eigen::VectorXd& xj = samples->at(it.first)->state[0];

    // check the neighbors
    const double diffnorm = (x-xj).norm();
    EXPECT_TRUE(diffnorm<r+1.0e-10);
    EXPECT_NEAR(diffnorm*diffnorm, it.second, 1.0e-10);
  }
}

TEST_F(GraphLaplacianTests, FindNearestNeighbors_NumNeighbors) {
  // create the graph laplacian from samples
  auto samples = CreateFromSamples();

  // build the kd-tree based on the samples
  laplacian->BuildKDTree();

  // choose a new random point from the distribution
  const Eigen::VectorXd x = rv->Sample();
  const size_t k = 15; // the number of nearest neighbors

  // find all of the points within a choosen radius
  std::vector<std::pair<std::size_t, double> > neighbors;
  laplacian->FindNeighbors(x, k, neighbors);

  // make sure we found the correct number of nearest neighbors
  EXPECT_EQ(neighbors.size(), k);

  // make sure the index of the neighbors is valid and they are within the correct radius
  for( const auto& it : neighbors ) {
    EXPECT_TRUE(it.first<n);

    // make sure that we have found the correct neighbor
    const Eigen::VectorXd& xj = samples->at(it.first)->state[0];

    // check the neighbors
    const double diffnorm = (x-xj).norm();
    EXPECT_NEAR(diffnorm*diffnorm, it.second, 1.0e-10);
  }
}
