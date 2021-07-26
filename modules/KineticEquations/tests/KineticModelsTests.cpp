#include <gtest/gtest.h>

#include "spipack/KineticEquations/KineticModels.hpp"

using namespace spi::KineticEquations;

class CustomKineticModelsTest : public KineticModels {
public:

  CustomKineticModelsTest() : KineticModels() {}

  inline static double InitialMassDensity(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x) { return x(0); }

  inline static Eigen::Matrix<double, MacroscaleInformation::dim, 1> InitialExpectedVelocity(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x) { return Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Ones(MacroscaleInformation::dim, 1); }

  inline static double InitialExpectedEnergy(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x) { return 2.0; }

private:

};

TEST(KineticModelsTests, StaticDefaultInitialConditions) {
  // chose a point in the domain
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> x = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Constant(MacroscaleInformation::dim, 0.5);

  // check initial mass density with gradient
  EXPECT_DOUBLE_EQ(KineticModels::InitialMassDensity(x), 1.0);

  // check initial expected velocity and its divergence
  const Eigen::Matrix<double, MacroscaleInformation::dim, 1> vel = KineticModels::InitialExpectedVelocity(x);
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) { EXPECT_DOUBLE_EQ(vel(i), 0.0); }

  // check initial expected energy
  EXPECT_DOUBLE_EQ(KineticModels::InitialExpectedEnergy(x), 1.0);
}

TEST(KineticModelsTests, DefaultInitialConditions) {
  // construct the kinetic models
  auto kinetic = std::make_shared<KineticModels>();
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

TEST(KineticModelsTests, CustomInitialConditions) {
  // construct the kinetic models
  auto kinetic = std::make_shared<CustomKineticModelsTest>();
  assert(kinetic);

  // chose a point in the domain
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> x = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Constant(MacroscaleInformation::dim, 0.5);

  // check initial mass density with gradient
  EXPECT_DOUBLE_EQ(kinetic->InitialMassDensity(x), x(0));

  // check initial expected velocity and its divergence
  const Eigen::Matrix<double, MacroscaleInformation::dim, 1> vel = kinetic->InitialExpectedVelocity(x);
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) { EXPECT_DOUBLE_EQ(vel(i), 1.0); }

  // check initial expected energy
  EXPECT_DOUBLE_EQ(kinetic->InitialExpectedEnergy(x), 2.0);
}
