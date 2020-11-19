#include <gtest/gtest.h>

#include "spipack/ConservationEquations/CompressibleEuler.hpp"

using namespace spi::ConservationEquations;

TEST(CompressibleEulerTests, Construct) {
  CompressibleEuler ceEquation;

  std::cout << "DONE" << std::endl;
}
