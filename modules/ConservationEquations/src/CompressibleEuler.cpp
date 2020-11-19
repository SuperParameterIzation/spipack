#include "spipack/ConservationEquations/CompressibleEuler.hpp"

#include <iostream>
#include <vector>

#include <exahype/main.h>

using namespace spi::ConservationEquations;

CompressibleEuler::CompressibleEuler() {
  char** st = new char*[2];
  st[0] = ""; st[1] = "/home/andy/Software/SuperParameterIzation.github.io/spipack/spipack/ConservationEquations/ExaHyPE/GenerateCompressibleEuler.exahype";

  exahype::main(2, st);

  delete[] st;


}
