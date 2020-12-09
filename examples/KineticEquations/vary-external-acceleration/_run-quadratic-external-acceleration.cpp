#include "QuadraticExternalAcceleration.hpp"

int main(int argc, char **argv) {
  // the output filename
  const std::string filename = "output/quadratic-external-acceleration/output";

  // create the conditional velocity distribution
  auto distribution = std::make_shared<QuadraticExternalAcceleration>(VaryExternalAcceleration::Options(filename));

  // run the micro-scale model
  double time = 0.0;
  const double delta = 0.1;
  for( std::size_t i=0; i<5; ++i ) {
    time += delta;
    distribution->Run(time, VaryExternalAcceleration::DefaultMacroscaleInformation(), i==0);
  }
}
