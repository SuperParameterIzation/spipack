#include <gtest/gtest.h>

#include "spipack/KineticEquations/MacroscaleInformation.hpp"

using namespace spi::KineticEquations;

TEST(MacroscaleInformationTests, DefaultSetup) {
  // the spatial dimension is 2
  EXPECT_EQ(MacroscaleInformation::dim, 2);

  // create the macro-scale information with defaults
  auto macroInfo = std::make_shared<MacroscaleInformation>();

  // the default mass density is 1
  EXPECT_DOUBLE_EQ(macroInfo->massDensity, 1.0);

  // the default log mass density gradient is 0
  EXPECT_EQ(macroInfo->logMassDensityGrad.size(), MacroscaleInformation::dim);
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
    EXPECT_DOUBLE_EQ(macroInfo->logMassDensityGrad(i), 0.0);
  }

  // the default velocity is 0
  EXPECT_EQ(macroInfo->velocity.size(), MacroscaleInformation::dim);
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
    EXPECT_DOUBLE_EQ(macroInfo->velocity(i), 0.0);
  }

  // test the coordinate system
  EXPECT_TRUE(macroInfo->coordinates==Coordinates::MACRO);
}

TEST(MacroscaleInformationTests, Setup) {
  // the spatial dimension is 2
  EXPECT_EQ(MacroscaleInformation::dim, 2);

  const double massDensity = 1.5;
  const Eigen::Matrix<double, MacroscaleInformation::dim, 1> logMassDensityGrad = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Random(MacroscaleInformation::dim, 1);
  const Eigen::Matrix<double, MacroscaleInformation::dim, 1> velocity = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Random(MacroscaleInformation::dim, 1);
  const double velocityDiv = 1.0;

  // create the macro-scale information
  auto macroInfo = std::make_shared<MacroscaleInformation>(massDensity, logMassDensityGrad, velocity, velocityDiv, Coordinates::MICRO);

  // check the mass density
  EXPECT_DOUBLE_EQ(macroInfo->massDensity, massDensity);

  // check the log mass density gradient
  EXPECT_EQ(macroInfo->logMassDensityGrad.size(), MacroscaleInformation::dim);
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
    EXPECT_DOUBLE_EQ(macroInfo->logMassDensityGrad(i), logMassDensityGrad(i));
  }

  // check the velocity
  EXPECT_EQ(macroInfo->velocity.size(), MacroscaleInformation::dim);
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
    EXPECT_DOUBLE_EQ(macroInfo->velocity(i), velocity(i));
  }

  // test the coordinate system
  EXPECT_TRUE(macroInfo->coordinates==Coordinates::MICRO);
}
