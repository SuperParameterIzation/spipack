#include <gtest/gtest.h>

#include <MUQ/Utilities/RandomGenerator.h>
#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/KineticEquations/KineticParticleModels.hpp"

#include "spipack/ConservationEquations/Boltzmann.hpp"

using namespace muq::Modeling;
using namespace spi::KineticEquations;
using namespace spi::ConservationEquations;

/*class Work {
public:
  Work(double aIn) : a(aIn){};

  double Evaluate(double x) {
    return a*x*x;
  };

  double AnotherEvaluate(double x) {
    return a*a*x+x*x-a;
  };

private:
  const double a;
};*/

TEST(ParallelBoltzmanTests, BoltzmannParticleMethods) {
  // set up the communicator
  auto comm = std::make_shared<parcer::Communicator>();
  const int rank = comm->GetRank();

  // the number of particles per micro-scale models
  const std::size_t n = 1500;

  // the number of timesteps per micro-scale model run
  const std::size_t numTimesteps = 5;

  // number of micro-scale models
  std::size_t gridSize = 13;

  // options for the nearest neighbor search
  YAML::Node nnOptions;
  nnOptions["NumSamples"] = n;
  nnOptions["Stride"] = n/5;
  nnOptions["NumThreads"] = 1;

  // options for the bandwidth parameter optimization within the Kolmogorov problem
  YAML::Node kolParaOptimization;
  kolParaOptimization["SparsityTolerance"] = 1.0e-2;
  kolParaOptimization["FTol.AbsoluteTolerance"] = 1.0e-2;
  kolParaOptimization["FTol.RelativeTolerance"] = 1.0e-2;
  kolParaOptimization["XTol.AbsoluteTolerance"] = 1.0e-2;
  kolParaOptimization["XTol.RelativeTolerance"] = 1.0e-2;
  kolParaOptimization["MaxEvaluations"] = 1000;
  kolParaOptimization["Algorithm"] = "LBFGS";

  // options for the Kolmogorov operator
  YAML::Node kolOptions;
  kolOptions["EigensolverTolerance"] = 1.0e-5;
  kolOptions["TruncationTolerance"] = -std::log(1.0e-4);
  kolOptions["NumNearestNeighbors"] = std::min((std::size_t)25, n/2);
  kolOptions["NumEigenvalues"] = (std::size_t)(5*log((double)n));
  kolOptions["BandwidthCostOptimization"] = kolParaOptimization;

  // options for the conditional velocity distribution
  YAML::Node conditionalOptions;
  conditionalOptions["NearestNeighbors"] = nnOptions;
  conditionalOptions["NumTimesteps"] = numTimesteps;
  conditionalOptions["KolmogorovOptions"] = kolOptions;
  conditionalOptions["AccelerationNoiseScale"] = 0.0;
  conditionalOptions["NondimensionalParameter"] = 1.0e-2;

  // set the options for the kinetic models
  YAML::Node options;
  options["GridSize"] = gridSize;
  options["ConditionalVelocityDistribution"] = conditionalOptions;

  // create a standard Gaussian random variable
  auto rv = std::make_shared<Gaussian>(MacroscaleInformation::dim)->AsVariable();

  // construct the micro-scale models
  auto kinetic = KineticParticleModels::Construct(rv, options, comm);
  assert(kinetic);

  if( rank==0 ) {
    // the kinetic models are not yet registered
    EXPECT_FALSE(KineticModels::ExaHyPEKineticModels());

    // create the solver
    Boltzmann solver(kinetic, options);

    // make sure the kinetic models are registered
    EXPECT_TRUE(KineticModels::ExaHyPEKineticModels());

    solver.Run();
  }

  std::cout << "rank " << rank << " finished" << std::endl;


  /*// Set up the work to perform (this is done on all processes)
  std::shared_ptr<Work> work = std::make_shared<Work>(2.0);

  {
    // Set up the Queue.  This will "capture" all processes but the master (rank=0)
    auto comm = std::make_shared<parcer::Communicator>();

    // Create the queue.  Template parameters are <InputType, OutputType, WorkType>
    parcer::Queue<double, double, Work> queue(work, comm);

    //Use the master process to submit jobs.
    if(rank==0){
      unsigned numEvals = 100;
      std::vector<unsigned> ptIds(numEvals);
      std::vector<double> evalPts(numEvals);

      // Submit all the work to the queue
      for(unsigned i=0; i<numEvals; ++i){
        evalPts.at(i) = double(i)/double(numEvals-1);
        ptIds.at(i) = queue.SubmitWork(evalPts.at(i));
      }

      // Retrieve all the completed evaluations. Each call will hang until that evaluation is available
      for(unsigned i=0; i<numEvals; ++i){
        double eval = queue.GetResult(ptIds.at(i));
        EXPECT_DOUBLE_EQ(work->Evaluate(evalPts.at(i)), eval);
      }
    }
  }
  */
}
