#include "spipack/ConservationEquations/Boltzmann.hpp"

#include <iostream>

#include <exahype/main.h>

using namespace spi::ConservationEquations;

Boltzmann::Boltzmann() {
  std::cout << "HERE" << std::endl;

  char** st = new char*[2];
  st[0] = "";
  st[1] = "/home/andy/Software/SuperParameterIzation.github.io/spipack/spipack/ConservationEquations/ExaHyPE/BoltzmannEuler.exahype";

  exahype::main(2, st);

  delete[] st[0]; delete[] st[1];
  delete[] st;

}
