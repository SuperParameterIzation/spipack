#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/KineticEquations/KineticParticleModels.hpp"

using namespace muq::Modeling;
using namespace spi::KineticEquations;

TEST(KineticParticleModelsTests, DefaultInitialConditions) {
  // construct the kinetic models
  auto kinetic = std::make_shared<KineticParticleModels>();
  assert(kinetic);

  // chose a point in the domain
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> x = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Constant(MacroscaleInformation::dim, 0.5);

  // check initial mass density with gradient
  EXPECT_DOUBLE_EQ(kinetic->InitialMassDensity(x), 1.0);

  // check initial expected velocity and its divergence
  const Eigen::Matrix<double, MacroscaleInformation::dim, 1> vel = kinetic->InitialExpectedVelocity(x);
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) { EXPECT_DOUBLE_EQ(vel(i), 0.0); }

  // check initial expected energy
  EXPECT_DOUBLE_EQ(kinetic->InitialExpectedEnergy(x), 1.0);
}

TEST(KineticParticleModelsTests, ModelGridConstruction) {
  // the number of particles per micro-scale models
  const std::size_t n = 1000;

  // the number of timesteps per micro-scale model run
  const std::size_t numTimesteps = 5;

  // number of micro-scale models
  std::size_t gridSize = 10;

  // options for the nearest neighbor search
  YAML::Node nnOptions;
  nnOptions["NumSamples"] = n;
  nnOptions["Stride"] = n/5;
  nnOptions["NumThreads"] = 1;

  // options for the Kolmogorov operator
  YAML::Node kolOptions;
  kolOptions["EigensolverTolerance"] = 1.0e-5;
  kolOptions["TruncationTolerance"] = -std::log(1.0e-4);
  kolOptions["NumNearestNeighbors"] = std::min((std::size_t)25, n/2);
  kolOptions["NumEigenvalues"] = (std::size_t)(5*log((double)n));

  // options for the conditional velocity distribution
  YAML::Node conditionalOptions;
  conditionalOptions["NearestNeighbors"] = nnOptions;
  conditionalOptions["NumTimesteps"] = numTimesteps;
  conditionalOptions["KolmogorovOptions"] = kolOptions;
  conditionalOptions["AccelerationNoiseScale"] = 0.0;
  conditionalOptions["NondimensionalParameter"] = 1.0e-6;

  // set the options for the kinetic models
  YAML::Node options;
  options["GridSize"] = gridSize;
  options["ConditionalVelocityDistribution"] = conditionalOptions;

  // create a random variable to sample as initial conditions
  auto rv = std::make_shared<Gaussian>(MacroscaleInformation::dim)->AsVariable();

  // create the kinetic models
  auto kinetic = KineticParticleModels::Construct(rv, options);

  // check the grid size
  EXPECT_EQ(kinetic->GridSize(), gridSize);
  EXPECT_EQ(kinetic->NumMicroscaleModels(), gridSize*gridSize);

  // check each micro-scale model
  for( std::size_t i=0; i<kinetic->NumMicroscaleModels(); ++i ) {
    EXPECT_TRUE(kinetic->MicroscaleModel(i));
    EXPECT_EQ(kinetic->MicroscaleModel(i)->NumSamples(), n);
  }

  // check each micro-scale model (iterators)
  for( auto it=kinetic->ModelIteratorBegin(); it!=kinetic->ModelIteratorEnd(); ++it ) {
    EXPECT_TRUE(*it);
    EXPECT_EQ((*it)->NumSamples(), n);
  }
}
