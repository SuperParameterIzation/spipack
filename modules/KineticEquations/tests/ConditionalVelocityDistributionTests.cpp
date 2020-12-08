#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/KineticEquations/ConditionalVelocityDistribution.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::KineticEquations;

class ConditionalVelocityDistributionTests : public::testing::Test {
protected:

  /// Set up information to test the sample representation
  virtual inline void SetUp() override {
    // create a standard Gaussian random variable
    rv = std::make_shared<Gaussian>(dim)->AsVariable();

    // options for the nearest neighbor search
    YAML::Node nnOptions;
    nnOptions["NumSamples"] = n;
    nnOptions["Stride"] = n/5;
    nnOptions["NumThreads"] = omp_get_max_threads();

    // set the options for the conditional velocity distribution
    options["NearestNeighbors"] = nnOptions;
    options["NumTimesteps"] = numTimesteps;

    // construct the initial macro-scale information
    macroInfo = std::make_shared<ConditionalVelocityDistribution::MacroscaleInformation>(massDensity, expectedVel, expectedVelDiv, logMassDensityGrad);
  }

  /// Make sure everything is what we expect
  virtual inline void TearDown() override {
    // check the parameters
    EXPECT_EQ(distribution->NumSamples(), n);
    EXPECT_DOUBLE_EQ(distribution->TimestepParameter(), 0.5);
    EXPECT_EQ(distribution->NumTimesteps(), numTimesteps);
    EXPECT_DOUBLE_EQ(distribution->ExternalAccelerationRescaling(), 1.0);

    // check the macro-scale information
    auto macro = distribution->MacroscaleInfo();
    EXPECT_DOUBLE_EQ(macro->velocityDivergence, expectedVelDiv);
  }

  /// Create the conditional velocity distribution from samples
  /**
    \return The sample collection used to create the sample representation
  */
  inline std::shared_ptr<SampleCollection> CreateFromSamples() {
    // add random samples into a sample collection
    auto samples = std::make_shared<SampleCollection>();
    for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

    // create the sample representation
    distribution = std::make_shared<ConditionalVelocityDistribution>(macroLoc, samples, macroInfo, options);

    // return the samples
    return samples;
  }

  /// The dimension of state spaces
  const unsigned int dim = 4;

  /// The number of samples
  const std::size_t n = 1000;

  /// The number of timesteps
  const std::size_t numTimesteps = 100;

  /// The mass density
  double massDensity = 1.0;

  /// The expected velocity
  Eigen::VectorXd expectedVel = Eigen::VectorXd::Zero(dim);

  /// The expected velocity divergence
  double expectedVelDiv = 0.0;

  /// The gradient of the log-mass density
  Eigen::VectorXd logMassDensityGrad = Eigen::VectorXd::Zero(dim);

  /// The macroscale location
  Eigen::VectorXd macroLoc = Eigen::VectorXd::Zero(dim);

  /// The random variable that we sample for the initial conditions
  std::shared_ptr<RandomVariable> rv;

  /// The options for the graph Laplacian
  YAML::Node options;

  /// The initial macro-scale information
  std::shared_ptr<ConditionalVelocityDistribution::MacroscaleInformation> macroInfo;

  /// The conditional velocity distribution---use a pointer here so we can initalize it as null
  std::shared_ptr<ConditionalVelocityDistribution> distribution;
};

TEST_F(ConditionalVelocityDistributionTests, RandomVariableConstruction) {
  // create the conditional velocity distribution
  distribution = std::make_shared<ConditionalVelocityDistribution>(macroLoc, rv, macroInfo, options);

  // check the macro-scale information
  auto macro = distribution->MacroscaleInfo();
  EXPECT_DOUBLE_EQ(macro->massDensity, macroInfo->massDensity);
  EXPECT_NEAR((macro->velocity-expectedVel).norm(), 0.0, 1.0e-10);
  EXPECT_DOUBLE_EQ(macro->velocityDivergence, macroInfo->velocityDivergence);
  EXPECT_NEAR((macro->logMassDensityGrad-logMassDensityGrad).norm(), 0.0, 1.0e-10);

  // check the normalizing constant
  EXPECT_DOUBLE_EQ(distribution->NormalizingConstant(), 1.0);
}

TEST_F(ConditionalVelocityDistributionTests, SampleCollectionConstruction) {
  // create the Kolmogorov operator from samples
  auto samples = CreateFromSamples();

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-distribution->Point(i)).norm(), 0.0, 1.0e-10);
  }

  // check the macro-scale information
  auto macro = distribution->MacroscaleInfo();
  EXPECT_DOUBLE_EQ(macro->massDensity, macroInfo->massDensity);
  EXPECT_NEAR((macro->velocity-expectedVel).norm(), 0.0, 1.0e-10);
  EXPECT_DOUBLE_EQ(macro->velocityDivergence, macroInfo->velocityDivergence);
  EXPECT_NEAR((macro->logMassDensityGrad-logMassDensityGrad).norm(), 0.0, 1.0e-10);

  // check the normalizing constant
  EXPECT_DOUBLE_EQ(distribution->NormalizingConstant(), 1.0);
}

TEST_F(ConditionalVelocityDistributionTests, NearestNeighborsConstruction) {
  // add random samples into a sample collection
  auto samples = std::make_shared<SampleCollection>();
  for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

  auto neighbors = std::make_shared<NearestNeighbors>(samples, options["NearestNeighbors"]);

  // create the sample representation
  distribution = std::make_shared<ConditionalVelocityDistribution>(macroLoc, neighbors, macroInfo, options);

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-distribution->Point(i)).norm(), 0.0, 1.0e-10);
  }

  // check the macro-scale information
  auto macro = distribution->MacroscaleInfo();
  EXPECT_DOUBLE_EQ(macro->massDensity, macroInfo->massDensity);
  EXPECT_NEAR((macro->velocity-expectedVel).norm(), 0.0, 1.0e-10);
  EXPECT_DOUBLE_EQ(macro->velocityDivergence, macroInfo->velocityDivergence);
  EXPECT_NEAR((macro->logMassDensityGrad-logMassDensityGrad).norm(), 0.0, 1.0e-10);

  // check the normalizing constant
  EXPECT_DOUBLE_EQ(distribution->NormalizingConstant(), 1.0);
}

TEST_F(ConditionalVelocityDistributionTests, RunTest) {
  // create the conditional velocity distribution
  distribution = std::make_shared<ConditionalVelocityDistribution>(macroLoc, rv, macroInfo, options);

  // the current time is (by default) set to zero
  EXPECT_NEAR(distribution->CurrentTime(), 0.0, 1.0e-12);

  // construct the macro-scale information
  massDensity = 1.25;
  expectedVel = Eigen::VectorXd::Zero(dim);
  expectedVel(0) = 1.0;
  expectedVelDiv = 1.0;
  logMassDensityGrad = Eigen::VectorXd::Zero(dim);
  logMassDensityGrad(0) = 1.0;
  auto nextMacroInfo = std::make_shared<ConditionalVelocityDistribution::MacroscaleInformation>(massDensity, expectedVel, expectedVelDiv, logMassDensityGrad);

  // run the micro-scale model
  const double nextTime = 0.1;
  distribution->Run(nextTime, nextMacroInfo);
  EXPECT_NEAR(distribution->CurrentTime(), nextTime, 1.0e-12);
}
