#include "spipack/Tools/Kernels/HatKernel.hpp"

using namespace spi::Tools;

HatKernel::HatKernel() : CompactKernel() {}

double HatKernel::EvaluateCompactKernel(double const theta) const {
  // we have already checked that theta<=1---this function implements the support of the kernel. Therefore, we just need to return 1
  return 1.0;
}
