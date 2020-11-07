#include "spipack/Tools/Kernels/BumpKernel.hpp"

using namespace spi::Tools;

SPIPACK_REGISTER_KERNEL(BumpKernel)
SPIPACK_REGISTER_ISOTROPIC_KERNEL(BumpKernel)
SPIPACK_REGISTER_COMPACT_KERNEL(BumpKernel)

BumpKernel::BumpKernel(YAML::Node const& options) :
  CompactKernel(options),
  mag(options["Magnitude"].as<double>(1.0)),
  scale(options["Scale"].as<double>(1.0)),
  expon(options["Exponent"].as<double>(1.0))
{}

double BumpKernel::EvaluateCompactKernel(double const theta) const {
  // if theta is near 1, return zero to avoid a seg fault
  if( std::abs(theta-1.0)<1.0e-10 ) { return 0.0; }

  return mag*std::exp(scale*(1.0-1.0/(1.0-std::pow(theta, expon))));
}
