#include "spipack/Tools/Kernels/HatKernel.hpp"

#include <tgmath.h>

using namespace spi::Tools;

SPIPACK_REGISTER_KERNEL(HatKernel)
SPIPACK_REGISTER_ISOTROPIC_KERNEL(HatKernel)
SPIPACK_REGISTER_COMPACT_KERNEL(HatKernel)

HatKernel::HatKernel(YAML::Node const& options) :
  CompactKernel(options),
  mag(options["Magnitude"].as<double>(1.0))
{}

double HatKernel::EvaluateCompactKernelImpl(double const theta) const { return mag; }

double HatKernel::IsotropicKernelDerivative(double const theta) const { return 0.0; }

double HatKernel::IsotropicKernelSecondDerivative(double const theta) const { return 0.0; }

double HatKernel::Magnitude() const { return mag; }

double HatKernel::Integrate(std::size_t const dim, std::size_t const n) const {
  if( n==0 ) {
    return mag*std::pow(M_PI, dim/2.0)/std::tgamma(dim/2.0+1.0)/2.0;
  }

  return NumericallyIntegrate(dim, n);
}
