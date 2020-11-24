#include "spipack/ConservationEquations/Boltzmann.hpp"

#include "spipack/SPIPACKDirectories.hpp"





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

Boltzmann::Boltzmann() :
exahype::runners::Runner(parser, cmdlineargs)
{
  // turn on/off ExaHyPE logging
  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  tarch::logging::LogFilterFileReader::parsePlainTextFile(spi::SPIPACK_INSTALL_DIR+"/ExaHyPE/exahype.log-filter");
}

void Boltzmann::ConstructBoltzmannParser(std::stringstream& specfile) {
  // the file name of the ExaHyPE configuration
  const std::string exahypeConfig = spi::SPIPACK_INSTALL_DIR+"/ExaHyPE/GenerateBoltzmann.exahype";

  // read the input file into a parser object
  specfile.str(kernels::readSpecificationFileToJSON(exahypeConfig));
  parser.readFile(specfile, exahypeConfig);

  // check to make sure it read the correct input file
  if( !parser.isValid() ) {
    std::cerr << "ERROR: Invalid config file when constructing Boltzmann solver." << std::endl;
    assert(false);
  }
}

void Boltzmann::Run() {
  // setup environment
  peano::fillLookupTables();

  // parse specification file
  std::stringstream specfile;
  ConstructBoltzmannParser(specfile);

  // init solver registries
  kernels::registerSolvers(parser);

  // configure the logging
  tarch::logging::CommandLineLogger::getInstance().clearFilterList();

  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      false,  // logMachineName
      true,   // logMessageType
      false,   // logTrace
      parser.getLogFileName() );

  // turn on/off ExaHyPE logging
  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  tarch::logging::LogFilterFileReader::parsePlainTextFile(spi::SPIPACK_INSTALL_DIR+"/ExaHyPE/exahype.log-filter");

  assert(parser.isValid());
  int result = 0;

  initOptimisations();
  initHeaps();

  auto* repository = createRepository();

  // must come after repository creation
  initSolvers();

  if( parser.getMeasureCellProcessingTimes() ) { measureCellProcessingTimes(); }

  initDistributedMemoryConfiguration();
  initSharedMemoryConfiguration();
  initDataCompression();
  initHPCEnvironment();

  Run(*repository);

  if ( _parser.isValid() )
  shutdownSharedMemoryConfiguration();
  if ( _parser.isValid() )
  shutdownDistributedMemoryConfiguration();

  shutdownHeaps();

  delete repository;

  peano::shutdownParallelEnvironment();
  peano::shutdownSharedMemoryEnvironment();
  peano::releaseCachedData();

  kernels::finalise();
}

void Boltzmann::Run(exahype::repositories::Repository& repository) {
  peano::utils::UserInterface::writeHeader();

  if( !exahype::solvers::RegisteredSolvers.empty() ) {
    initialiseMesh(repository);

    bool communicatePeanoVertices =
        !exahype::solvers::Solver::DisablePeanoNeighbourExchangeInTimeSteps;

    if ( !_parser.getProfileEmptyAdapter() ) {
      repository.switchToInitialPrediction();
      repository.iterate( exahype::solvers::Solver::PredictionSweeps, communicatePeanoVertices );
      logInfo("runAsMaster(...)","computed first predictor (number of predictor sweeps: "<<exahype::solvers::Solver::PredictionSweeps<<")");
    }

    // configure time stepping loop
    double simulationEndTime = std::numeric_limits<double>::infinity();
    int simulationTimeSteps  = std::numeric_limits<int>::max();
    if (_parser.foundSimulationEndTime()) {
      simulationEndTime = _parser.getSimulationEndTime();
    } else {
      simulationTimeSteps = _parser.getSimulationTimeSteps();
    }
    const bool profileEmptyAdapter = _parser.getProfileEmptyAdapter();
    const exahype::parser::Parser::ProfilingTarget profilingTarget = _parser.getProfilingTarget();
    if ( profilingTarget==exahype::parser::Parser::ProfilingTarget::WholeCode ) {
      printTimeStepInfo(-1,repository);
      validateInitialSolverTimeStepData(exahype::solvers::Solver::FuseAllADERDGPhases);
    }
    const bool skipReductionInBatchedTimeSteps          = _parser.getSkipReductionInBatchedTimeSteps();
    const bool fuseMostADERDGPhases                     = _parser.getFuseMostAlgorithmicSteps();
    const bool fuseMostADERDGPhasesDoRiemannSolvesTwice = _parser.getFuseMostAlgorithmicStepsDoRiemannSolvesTwice();

    // run time stepping loop
    int timeStep = 0;
    while (
        tarch::la::greater(exahype::solvers::Solver::getMinTimeStepSizeOfAllSolvers(), 0.0) &&
        exahype::solvers::Solver::getMinTimeStampOfAllSolvers() < simulationEndTime         &&
        timeStep < simulationTimeSteps
    ) {
      UpdateSmallScale();

      bool plot = exahype::plotters::checkWhetherPlotterBecomesActive(
          exahype::solvers::Solver::getMinTimeStampOfAllSolvers()); // has no side effects

      preProcessTimeStepInSharedMemoryEnvironment();

      // fused time stepping
      int numberOfStepsToRun = 1;
      if ( profilingTarget==exahype::parser::Parser::ProfilingTarget::WholeCode && exahype::solvers::Solver::FuseAllADERDGPhases  ) {
        if (plot) {
          numberOfStepsToRun = 0;
        }
        else if (
            skipReductionInBatchedTimeSteps &&
            exahype::solvers::Solver::allSolversUseTimeSteppingScheme(exahype::solvers::Solver::TimeStepping::GlobalFixed) &&
            repository.getState().getVerticalExchangeOfSolverDataRequired()==false // known after mesh update
        ) {
          numberOfStepsToRun = determineNumberOfBatchedTimeSteps(timeStep);
        }

        runTimeStepsWithFusedAlgorithmicSteps(repository,numberOfStepsToRun);
        printTimeStepInfo(numberOfStepsToRun,repository);
      // straightforward realisation
    } else if ( profilingTarget==exahype::parser::Parser::ProfilingTarget::WholeCode ) {
        if ( fuseMostADERDGPhases ) {
          runOneTimeStepWithTwoSeparateAlgorithmicSteps(repository, plot);
        } else if ( fuseMostADERDGPhasesDoRiemannSolvesTwice ) {
          runOneTimeStepWithTwoSeparateAlgorithmicStepsDoRiemannSolvesTwice(repository, plot);
        } else {
          runOneTimeStepWithThreeSeparateAlgorithmicSteps(repository, plot);
        }
      }
      // profiling of isolated adapters
      else if ( profilingTarget==exahype::parser::Parser::ProfilingTarget::Prediction ) {
        logInfo("runAsMaster(...)","step " << timeStep << "\t\trun "<<exahype::solvers::Solver::PredictionSweeps<<" iteration with PredictionRerun adapter");
        printGridStatistics(repository);
        runPredictionInIsolation(repository);
      } else if ( profilingTarget==exahype::parser::Parser::ProfilingTarget::NeigbhourMerge ) {
        logInfo("runAsMaster(...)","step " << timeStep << "\t\trun one iteration with MergeNeighours adapter");
        printGridStatistics(repository);
        repository.switchToMergeNeighbours();
        repository.iterate(1,false);
      } else if ( profilingTarget==exahype::parser::Parser::ProfilingTarget::Update ) {
        logInfo("runAsMaster(...)","step " << timeStep << "\t\trun one iteration with UpdateAndReduce adapter");
        printGridStatistics(repository);
        repository.switchToUpdateAndReduce();
        repository.iterate(1,false);
      } else if ( profileEmptyAdapter ) {
        logInfo("runAsMaster(...)","step " << timeStep << "\t\trun one iteration with Empty adapter");
        printGridStatistics(repository);
        repository.switchToEmpty();
        repository.iterate(1,false);
      }

      postProcessTimeStepInSharedMemoryEnvironment();

      #if !defined(Parallel)
      logInfo("runAsMaster(...)", "memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
      #else
      if (tarch::parallel::Node::getInstance().getNumberOfNodes()==1) {
        logInfo("runAsMaster(...)", "memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
      }
      #endif

      timeStep += numberOfStepsToRun==0 ? 1 : numberOfStepsToRun;

      //std::cout << "BREAKING TIME LOOP!" << std::endl;
      //break;
    }
  }

  std::cout << "RAN SIMULATION!" << std::endl;

}

void Boltzmann::UpdateSmallScale() {
  std::cout << "time: " << exahype::solvers::Solver::getMinTimeStampOfAllSolvers() << std::endl;
  std::cout << "\tupdate small scale" << std::endl;

  //std::cout << offsetOfPatch[0] << std::endl;
}
