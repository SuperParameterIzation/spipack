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
    rv = std::make_shared<Gaussian>(MacroscaleInformation::dim)->AsVariable();

    velocity = Eigen::VectorXd::Zero(MacroscaleInformation::dim);
    logMassDensityGrad = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim);
    logMassDensityGrad(0) = 1.0;

    velocityDiv = 0.1;

    // options for the nearest neighbor search
    YAML::Node nnOptions;
    nnOptions["NumSamples"] = n;
    nnOptions["Stride"] = n/5;
    nnOptions["NumThreads"] = omp_get_max_threads();

    // options for the Kolmogorov operator
    static YAML::Node kolOptions;
    kolOptions["EigensolverTolerance"] = 1.0e-5;
    kolOptions["TruncationTolerance"] = -std::log(1.0e-4);
    kolOptions["NumNearestNeighbors"] = std::min((std::size_t)50, n/2);
    kolOptions["NumEigenvalues"] = (std::size_t)(5*log((double)n));

    // set the options for the conditional velocity distribution
    options["NearestNeighbors"] = nnOptions;
    options["KolmogorovOptions"] = kolOptions;
    options["NumTimesteps"] = numTimesteps;
    options["Verbosity"] = 2;

    // construct the initial macro-scale information using default values
    macroInfo = std::make_shared<MacroscaleInformation>(massDensity, logMassDensityGrad, velocity, velocityDiv, Coordinates::MACRO);
  }

  /// Make sure everything is what we expect
  virtual inline void TearDown() override {
    // check the parameters
    EXPECT_EQ(distribution->NumSamples(), n); // number of samples used to approximate the conditional distribution

    // check the macro-scale information
    auto macro = distribution->MacroscaleInfo();
    EXPECT_DOUBLE_EQ(macro->massDensity, macroInfo->massDensity);
    for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
      EXPECT_DOUBLE_EQ(macro->velocity(i), velocity(i));
      EXPECT_DOUBLE_EQ(macro->logMassDensityGrad(i), logMassDensityGrad(i));
    }

    // the sample mean should match the macro-scale velocity---the samples are shifted in the constructor so this is exact
    const Eigen::Matrix<double, MacroscaleInformation::dim, 1> sampleMeanVel = distribution->SampleMean();
    EXPECT_EQ(sampleMeanVel.size(), MacroscaleInformation::dim);
    for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
      EXPECT_NEAR(sampleMeanVel(i), velocity(i), 1.0e-10);
    }

    // check the macro-scale location
    const Eigen::Matrix<double, MacroscaleInformation::dim, 1> loc = distribution->MacroscaleLocation();
    EXPECT_EQ(loc.size(), MacroscaleInformation::dim);
    for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
      EXPECT_DOUBLE_EQ(loc(i), macroLoc(i));
    }
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

  /// The number of samples
  const std::size_t n = 1000;

  /// The number of timesteps
  const std::size_t numTimesteps = 5;

  /// The external acceleration rescaling
  const double alpha = 1.5;

  /// The mass density
  double massDensity = 1.25;

  /// The expected velocity
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> velocity = Eigen::VectorXd::Random(MacroscaleInformation::dim);

  /// The divergence of the expected velocity
  double velocityDiv = 0.1;

  /// The gradient of the log-mass density
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> logMassDensityGrad = Eigen::VectorXd::Random(MacroscaleInformation::dim);

  /// The macroscale location
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> macroLoc = Eigen::VectorXd::Random(MacroscaleInformation::dim);

  /// The random variable that we sample for the initial conditions
  std::shared_ptr<RandomVariable> rv;

  /// The options for the graph Laplacian
  YAML::Node options;

  /// The initial macro-scale information
  std::shared_ptr<MacroscaleInformation> macroInfo;

  /// The conditional velocity distribution---use a pointer here so we can initalize it as null
  std::shared_ptr<ConditionalVelocityDistribution> distribution;
};

TEST_F(ConditionalVelocityDistributionTests, RandomVariableConstruction) {
  // create the conditional velocity distribution
  distribution = std::make_shared<ConditionalVelocityDistribution>(macroLoc, rv, macroInfo, options);

  // make sure the kolmogorov operator is intialized
  //auto kol = distribution->LaplaceOperator();
  //EXPECT_TRUE(kol);
}

TEST_F(ConditionalVelocityDistributionTests, SampleCollectionConstruction) {
  // create the Kolmogorov operator from samples
  auto samples = CreateFromSamples();

  // check to make sure the samples match
  EXPECT_EQ(distribution->NumSamples(), n);
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-distribution->Point(i)).norm(), 0.0, 1.0e-10);
  }

  // make sure the kolmogorov operator is intialized
  /*auto kol = distribution->LaplaceOperator();
  EXPECT_TRUE(kol);
  EXPECT_EQ(kol->NumSamples(), n);
  for( std::size_t i=0; i<n; ++i ) {
    // after initialization, the kolmogorov samples are the exact same as the conditonal distribution samples (they should point to the same object)
    EXPECT_NEAR((samples->at(i)->state[0]-kol->Point(i)).norm(), 0.0, 1.0e-10);
  }*/
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

  // make sure the kolmogorov operator is intialized
  /*auto kol = distribution->LaplaceOperator();
  EXPECT_TRUE(kol);
  EXPECT_EQ(kol->NumSamples(), n);
  for( std::size_t i=0; i<n; ++i ) {
    // after initialization, the kolmogorov samples are the exact same as the conditonal distribution samples (they should point to the same object)
    EXPECT_NEAR((samples->at(i)->state[0]-kol->Point(i)).norm(), 0.0, 1.0e-10);
  }*/
}

TEST_F(ConditionalVelocityDistributionTests, RunTest) {
  options["OutputFilename"] = "ConditionalDistributionTest";
  options["NondimensionalParameter"] = 0.5;
  //options["NumTimesteps"] = 10;

  // create the conditional velocity distribution
  distribution = std::make_shared<ConditionalVelocityDistribution>(macroLoc, rv, macroInfo, options);

  // the current time is (by default) set to zero
  EXPECT_DOUBLE_EQ(distribution->CurrentTime(), 0.0);

  // construct the macro-scale information
  massDensity = 1.25;
  //logMassDensityGrad = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Random(MacroscaleInformation::dim);
  //velocity = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Random(MacroscaleInformation::dim);
  velocity = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim);
  velocity(0) = 1.0;

  // run the micro-scale model
  for( std::size_t t=1; t<75; ++t ) {
    const double nextTime = 0.01*t;
    velocity(0) = 0.1*t;

    auto nextMacroInfo = std::make_shared<MacroscaleInformation>(massDensity, logMassDensityGrad, velocity, velocityDiv, Coordinates::MACRO);

    distribution->Run(nextTime, nextMacroInfo);
    EXPECT_DOUBLE_EQ(distribution->CurrentTime(), nextTime);
  }

  /*// check the macro-scale information
  auto macro = distribution->MacroscaleInfo();
  EXPECT_DOUBLE_EQ(macro->massDensity, nextMacroInfo->massDensity);
  EXPECT_NEAR((macro->velocity-expectedVel).norm(), 0.0, 1.0e-10);
  EXPECT_NEAR((macro->logMassDensityGrad-logMassDensityGrad).norm(), 0.0, 1.0e-10);*/
}
