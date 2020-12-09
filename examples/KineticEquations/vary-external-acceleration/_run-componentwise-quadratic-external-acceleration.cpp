#include "ComponentwiseQuadraticExternalAcceleration.hpp"

int main(int argc, char **argv) {
  // the number of samples
  const std::size_t n = 1000;

  // the number of timesteps
  const std::size_t numTimesteps = 100;

  // the output filename
  std::string filename = "output/componentwise-quadratic-external-acceleration/output";

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

  // create the conditional velocity distribution
  auto distribution = std::make_shared<ComponentwiseQuadraticExternalAcceleration>(options);

  // run the micro-scale model
  double time = 0.0;
  const double delta = 0.1;
  for( std::size_t i=0; i<5; ++i ) {
    time += delta;
    distribution->Run(time, VaryExternalAcceleration::DefaultMacroscaleInformation());
  }
}
