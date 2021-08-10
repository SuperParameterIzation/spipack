#include "spipack/ConservationEquations/Boltzmann.hpp"

#include "spipack/ConservationEquations/ExaHyPE/Boltzmann/BoltzmannSolver.h"

#include "tarch/logging/Log.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "tarch/logging/CommandLineLogger.h"
#include "tarch/logging/LogFilterFileReader.h"
#include "tarch/parallel/Node.h"

#include "peano/peano.h"

#include "exahype/main.h"
#include "exahype/parser/Parser.h"
#include "exahype/Vertex.h"
#include "exahype/runners/Runner.h"
#include "exahype/repositories/Repository.h"
#include "exahype/plotters/Plotter.h"

#include "kernels/KernelCalls.h"

#include <vector>
#include <cstdlib> // getenv, exit
#include <iostream>
#include <cstdio>
#include "kernels/GaussLegendreBasis.h"

#include "peano/utils/UserInterface.h"

#include <iostream>

#include <exahype/main.h>

using namespace spi::ConservationEquations;

/*spiEX_Boltzmann::BoltzmannSolver* Boltzmann::BoltzmannSolver() {
  for( const auto& it : exahype::solvers::RegisteredSolvers ) {
    assert(it);
    auto bolt = dynamic_cast<spiEX_Boltzmann::BoltzmannSolver*>(it);
    if( bolt ) { return bolt; }
  }

  return nullptr;
}*/

void Boltzmann::Initialize() {
  // setup environment
  peano::fillLookupTables();
  int sharedMemorySetup = peano::initSharedMemoryEnvironment();
}

Boltzmann::Boltzmann(std::shared_ptr<spi::KineticEquations::KineticModels> const& microscaleModels, YAML::Node const& options) :
mode(SolverModesMap().at(options["SolverMode"].as<std::string>("CompressibleEulerLimit")))
{
  // setup environment
  peano::fillLookupTables();
  int sharedMemorySetup = peano::initSharedMemoryEnvironment();

  // turn on/off ExaHyPE logging
  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  tarch::logging::LogFilterFileReader::parsePlainTextFile(spi::SPIPACK_INSTALL_DIR+"/ExaHyPE/exahype.log-filter");

  // the file name of the ExaHyPE configuration
  std::vector<std::string> exahypeConfig(1, constructFile);
  std::stringstream specfile;
  specfile.str(kernels::readSpecificationFileToJSON(exahypeConfig[0]));
  exahype::parser::Parser parser;
  parser.readFile(specfile, exahypeConfig[0]);

  // init solver registries
  kernels::registerSolvers(parser);

  for( const auto& it : exahype::solvers::RegisteredSolvers ) {
    assert(it);
    auto bolt = dynamic_cast<spiEX_Boltzmann::BoltzmannSolver*>(it);
    if( bolt ) {
      bolt->kineticModels = microscaleModels;
      break;
    }
  }
}

Boltzmann::~Boltzmann() {
  // need to delete the particle methods so that the communicator releases the trapped processors (parallel only)
  for( const auto& it : exahype::solvers::RegisteredSolvers ) {
    assert(it);
    auto bolt = dynamic_cast<spiEX_Boltzmann::BoltzmannSolver*>(it);
    if( bolt ) {
      bolt->kineticModels = nullptr;
      break;
    }
  }
}

void Boltzmann::Run() {
  std::cout << "RUNNING THE BOLTZMANN EQUATION CODE" << std::endl;

  // the file name of the ExaHyPE configuration
  std::vector<std::string> exahypeConfig(1, constructFile);
  std::stringstream specfile;
  specfile.str(kernels::readSpecificationFileToJSON(exahypeConfig[0]));

  exahype::parser::Parser parser;
  parser.readFile(specfile, exahypeConfig[0]);
  assert(parser.isValid());

  exahype::runners::Runner runner(parser, exahypeConfig);
  int programExitCode = runner.run();

  peano::shutdownParallelEnvironment();
  peano::shutdownSharedMemoryEnvironment();
  peano::releaseCachedData();

  kernels::finalise();
}

Boltzmann::SolverMode Boltzmann::GetSolverMode() const { return mode; }

std::unordered_map<std::string, Boltzmann::SolverMode> Boltzmann::SolverModesMap() {
  // if we have already added the enums, just return the map
  static std::unordered_map<std::string, SolverMode> map;
  if( !map.empty() ) { return map; }

  // add the enum options to the map
  #define F(s) #s
  const std::string optionNames[] = { SPIPACK_BOLTZMANN_SOLVERMODES(F) };
  #undef F
  #define F(e) e
  const Boltzmann::SolverMode optionEnums[] = { SPIPACK_BOLTZMANN_SOLVERMODES(F) };
  #undef F

  for( std::size_t i=0; i<NumSolverModes; ++i ) { map[optionNames[i]] = optionEnums[i]; }

  return map;
}
