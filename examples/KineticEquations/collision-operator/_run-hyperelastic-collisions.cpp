#include "VaryCollisionEnergy.hpp"

int main(int argc, char **argv) {
  // the output filename
  const std::string filename = "output/hyperelastic-collisions/output";

  // create the conditional velocity distribution
  auto distribution = std::make_shared<VaryCollisionEnergy>(1.5, VaryCollisionEnergy::Options(filename));

  // run the micro-scale model
  double time = 0.0;
  const double delta = 0.05;
  for( std::size_t i=0; i<5; ++i ) {
    time += delta;
    distribution->Run(time, VaryCollisionEnergy::DefaultMacroscaleInformation(), i==0);
  }
}
