#include <gtest/gtest.h>

#include "spipack/ConservationEquations/Boltzmann.hpp"

using namespace spi::ConservationEquations;

TEST(BoltzmannTests, Construct) {
  Boltzmann boltzmannEquation;

  std::cout << "DONE" << std::endl;
}
