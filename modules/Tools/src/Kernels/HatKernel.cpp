#include "spipack/Tools/Kernels/HatKernel.hpp"

using namespace spi::Tools;

SPIPACK_REGISTER_KERNEL(HatKernel)
SPIPACK_REGISTER_ISOTROPIC_KERNEL(HatKernel)
SPIPACK_REGISTER_COMPACT_KERNEL(HatKernel)

HatKernel::HatKernel(YAML::Node const& options) :
  CompactKernel(options),
  mag(options["Magnitude"].as<double>(1.0))
{}

double HatKernel::EvaluateCompactKernel(double const theta) const {
  // we have already checked that theta<=1---this function implements the support of the kernel. Therefore, we just need to return the magnitude
  return mag;
}
