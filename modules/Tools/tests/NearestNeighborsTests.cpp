#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/Tools/NearestNeighbors.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;

class NearestNeighborsTests : public::testing::Test {
protected:
  /// Set up information to test the nearest neighbor construction
  virtual void SetUp() override {
    // create a standard Gaussian random variable
    rv = std::make_shared<Gaussian>(dim)->AsVariable();

    // set the options for the graph laplacian
    options["MaxLeaf"] = maxLeaf;
    options["NumSamples"] = n;
    options["Stride"] = 454;
  }

  /// Make sure everything is constructed correctly
  virtual void TearDown() override {
    EXPECT_TRUE(nn);
    EXPECT_TRUE(nn->Samples());

    // check the number of samples and the state dimension
    EXPECT_EQ(nn->NumSamples(), n);
    EXPECT_EQ(nn->StateDim(), dim);
  }

  /// Create the nearest neighbor searcher from samples
  /**
    \return The sample collection used to create the graph Laplacian
  */
  inline std::shared_ptr<SampleCollection> CreateFromSamples() {
    // add random samples into a sample collection
    auto samples = std::make_shared<SampleCollection>();
    assert(rv);
    for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

    // create the graph laplacian
    nn = std::make_shared<NearestNeighbors>(samples, options);

    // return the samples
    return samples;
  }

  /// The dimension of state spaces
  const std::size_t dim = 4;

  /// The number of samples
  const std::size_t n = 10000;

  /// The max leaf size for the kd tree
  const std::size_t maxLeaf = 15;

  /// The options for the graph Laplacian
  YAML::Node options;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The nearest neighbor searcher---use a pointer here so we can initalize it as null
  std::shared_ptr<NearestNeighbors> nn;
};

TEST_F(NearestNeighborsTests, RandomVariableConstruction) {
  nn = std::make_shared<NearestNeighbors>(rv, options);

  // construct the kd-trees
  nn->BuildKDTrees();

  // choose a random point
  const Eigen::VectorXd x = rv->Sample();

  // find the nearest neighbors to the point
  {
    const double radius = 0.1;
    std::vector<std::pair<std::size_t, double> > neighbors;
    nn->FindNeighbors(x, radius*radius, neighbors);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const double radius = 1.0;
    std::vector<std::pair<std::size_t, double> > neighbors;
    nn->FindNeighbors(x, radius*radius, neighbors, lag);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point
  {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    nn->FindNeighbors(x, k, neighbors);
    EXPECT_EQ(neighbors.size(), k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    nn->FindNeighbors(x, k, neighbors, lag);
    EXPECT_TRUE(neighbors.size()<=k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
    }
  }
}

TEST_F(NearestNeighborsTests, SampleCollectionConstruction) {
  auto samples = CreateFromSamples();
  EXPECT_TRUE(samples);

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-nn->Point(i)).norm(), 0.0, 1.0e-10);
  }

  // construct the kd-trees
  nn->BuildKDTrees();

  // choose a random point
  const Eigen::VectorXd x = rv->Sample();

  // find the nearest neighbors to the point
  {
    const double radius = 0.1;
    std::vector<std::pair<std::size_t, double> > neighbors;
    nn->FindNeighbors(x, radius*radius, neighbors);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const double radius = 1.0;
    std::vector<std::pair<std::size_t, double> > neighbors;
    nn->FindNeighbors(x, radius*radius, neighbors, lag);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
      EXPECT_TRUE(neigh.second<radius*radius);
    }
  }

  // find the nearest neighbors to the point
  {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    nn->FindNeighbors(x, k, neighbors);
    EXPECT_EQ(neighbors.size(), k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first<n);
    }
  }

  // find the nearest neighbors to the point with lag
  for( std::size_t lag=0; lag<n; lag+=10 ) {
    const std::size_t k = 10;
    std::vector<std::pair<std::size_t, double> > neighbors;
    nn->FindNeighbors(x, k, neighbors, lag);
    EXPECT_TRUE(neighbors.size()<=k);
    for( const auto& neigh : neighbors ) {
      EXPECT_TRUE(neigh.first>=lag);
      EXPECT_TRUE(neigh.first<n);
    }
  }
}

TEST_F(NearestNeighborsTests, SquaredBandwidth) {
  nn = std::make_shared<NearestNeighbors>(rv, options);

  // construct the kd-trees
  nn->BuildKDTrees();

  // compute the squared bandwidth
  const std::size_t k = 25; // number of nearest neighbors
  std::vector<std::vector<std::pair<std::size_t, double> > > neighbors;
  const Eigen::VectorXd squaredBandwidth = nn->SquaredBandwidth(k, neighbors);

  // loop through all of the points and compare the squared bandwidth to the expected squared bandwidth
  EXPECT_EQ(neighbors.size(), n);
  EXPECT_EQ(squaredBandwidth.size(), n);
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_EQ(neighbors[i].size(), k+1);

    double expected = 0.0;
    for( const auto& neigh : neighbors[i] ) { expected += neigh.second; }
    EXPECT_NEAR(squaredBandwidth(i), expected/k, 1.0e-10);
  }
}


TEST_F(NearestNeighborsTests, SquaredBandwidth_Multithreaded) {
  options["NumThreads"] = omp_get_max_threads();
  nn = std::make_shared<NearestNeighbors>(rv, options);
  EXPECT_EQ(nn->NumThreads(), omp_get_max_threads());

  // construct the kd-trees
  nn->BuildKDTrees();

  // compute the squared bandwidth
  const std::size_t k = 25; // number of nearest neighbors
  std::vector<std::vector<std::pair<std::size_t, double> > > neighbors;
  const Eigen::VectorXd squaredBandwidth = nn->SquaredBandwidth(k, neighbors);

  // loop through all of the points and compare the squared bandwidth to the expected squared bandwidth
  EXPECT_EQ(neighbors.size(), n);
  EXPECT_EQ(squaredBandwidth.size(), n);
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_EQ(neighbors[i].size(), k+1);

    double expected = 0.0;
    for( const auto& neigh : neighbors[i] ) { expected += neigh.second; }
    EXPECT_NEAR(squaredBandwidth(i), expected/k, 1.0e-10);
  }
}
