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

double BumpKernel::EvaluateCompactKernelImpl(double const theta) const {
  // if theta is near 1, return zero to avoid a seg fault
  if( std::abs(theta-1.0)<1.0e-10 ) { return 0.0; }

  return mag*std::exp(scale*(1.0-1.0/(1.0-std::pow(theta, expon))));
}

double BumpKernel::Magnitude() const { return mag; }

double BumpKernel::Scale() const { return scale; }

double BumpKernel::Exponent() const { return expon; }

double BumpKernel::IsotropicKernelDerivative(double const theta) const {
  // if theta is near 1, return zero to avoid a seg fault
  static const double thresh = 1.0-1.0e-10;
  if( theta>thresh ) { return 0.0; }

  const double powtheta = std::pow(theta, expon);
  const double _1mpowtheta = 1.0-powtheta;

  return -mag*scale*expon*powtheta/theta/(_1mpowtheta*_1mpowtheta)*std::exp(scale*(1.0-1.0/_1mpowtheta));
}

double BumpKernel::IsotropicKernelSecondDerivative(double const theta) const {
  // if theta is near 1, return zero to avoid a seg fault
  static const double thresh = 1.0-1.0e-10;
  if( theta>thresh ) { return 0.0; }

  const double powtheta = std::pow(theta, expon);
  const double _1mpowtheta = 1.0-powtheta;

  const double first = mag*std::exp(scale*(1.0-1.0/_1mpowtheta));
  const double firstderiv = -mag*scale*expon*powtheta/theta/(_1mpowtheta*_1mpowtheta)*std::exp(scale*(1.0-1.0/_1mpowtheta));

  const double top = -scale*expon*powtheta/theta;
  const double topderiv = -scale*expon*(expon-1.0)*powtheta/(theta*theta);
  const double bot = _1mpowtheta*_1mpowtheta;
  const double botderiv = -2.0*_1mpowtheta*expon*powtheta/theta;
  const double second = top/bot;
  const double seconderiv = (topderiv*bot-botderiv*top)/(bot*bot);

  return first*seconderiv + second*firstderiv;
}
