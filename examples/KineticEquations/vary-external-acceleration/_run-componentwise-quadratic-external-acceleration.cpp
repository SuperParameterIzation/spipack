#include "ComponentwiseQuadraticExternalAcceleration.hpp"

int main(int argc, char **argv) {
  // the output filename
  std::string filename = "output/componentwise-quadratic-external-acceleration/output";

  // create the conditional velocity distribution
  auto distribution = std::make_shared<ComponentwiseQuadraticExternalAcceleration>(VaryExternalAcceleration::Options(filename));

  // run the micro-scale model
  double time = 0.0;
  const double delta = 0.05;
  for( std::size_t i=0; i<5; ++i ) {
    time += delta;
    distribution->Run(time, VaryExternalAcceleration::DefaultMacroscaleInformation(), i==0);
  }
}
