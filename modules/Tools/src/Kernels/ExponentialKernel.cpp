#include "spipack/Tools/Kernels/ExponentialKernel.hpp"

using namespace spi::Tools;

SPIPACK_REGISTER_KERNEL(ExponentialKernel)
SPIPACK_REGISTER_ISOTROPIC_KERNEL(ExponentialKernel)

ExponentialKernel::ExponentialKernel(YAML::Node const& options) :
  IsotropicKernel(options),
  mag(options["Magnitude"].as<double>(1.0)),
  scale(options["Scale"].as<double>(1.0)),
  expon(options["Exponent"].as<double>(1.0))
{
  assert(mag>0.0);
  assert(scale>0.0);
  assert(expon>0.0);
}

double ExponentialKernel::EvaluateIsotropicKernel(double const theta) const { return mag*std::exp(-scale*std::pow(std::abs(theta), expon)); }

double ExponentialKernel::Magnitude() const { return mag; }

double ExponentialKernel::Scale() const { return scale; }

double ExponentialKernel::Exponent() const { return expon; }
