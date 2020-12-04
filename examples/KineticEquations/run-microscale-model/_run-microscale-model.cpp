#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/KineticEquations/ConditionalVelocityDistribution.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::KineticEquations;

int main(int argc, char **argv) {
  // the dimension of state spaces
  const unsigned int dim = 2;

  // create a standard Gaussian random variable
  auto rv = std::make_shared<Gaussian>(dim)->AsVariable();

  // the number of samples
  const std::size_t n = 1000;

  // the number of timesteps
  const std::size_t numTimesteps = 100;

  // the output filename
  std::string filename = "output/output";

  // options for the nearest neighbor search
  YAML::Node nnOptions;
  nnOptions["NumSamples"] = n;
  nnOptions["Stride"] = n/5;
  nnOptions["NumThreads"] = omp_get_max_threads();

  // set the options for the conditional velocity distribution
  YAML::Node options;
  options["NearestNeighbors"] = nnOptions;
  options["NumTimesteps"] = numTimesteps;
  options["OutputFilename"] = filename;

  // construct the initial macro-scale information
  auto macroInfo = std::make_shared<ConditionalVelocityDistribution::MacroscaleInformation>(0.0);

  // create the conditional velocity distribution
  auto distribution = std::make_shared<ConditionalVelocityDistribution>(rv, macroInfo, options);

  // construct the macro-scale information
  auto nextMacroInfo = std::make_shared<ConditionalVelocityDistribution::MacroscaleInformation>(1.0);

  // run the micro-scale model
  const double nextTime = 0.1;
  distribution->Run(nextTime, nextMacroInfo);
}
