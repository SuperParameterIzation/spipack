#include "spipack/Tools/Kernels/CompactKernel.hpp"

#include <iostream>

using namespace spi::Tools;

CompactKernel::CompactKernel(YAML::Node const& options) : IsotropicKernel(options) {}

double CompactKernel::EvaluateIsotropicKernel(double const theta) const {
  static const double buffer = 1.0+1.0e-10;
  return (theta<buffer? EvaluateCompactKernel(theta) : 0.0);
}
