#include "spipack/Tools/Kernels/IsotropicKernel.hpp"

#include <iostream>

using namespace spi::Tools;

IsotropicKernel::IsotropicKernel(YAML::Node const& options) : Kernel(options) {}

double IsotropicKernel::Evaluate(Eigen::Ref<Eigen::VectorXd> const& x1, Eigen::Ref<Eigen::VectorXd> const& x2) const {
  assert(x1.size()==x2.size());
  const Eigen::VectorXd diff = x1-x2;
  return EvaluateIsotropicKernel(diff.dot(diff));
}

double IsotropicKernel::operator()(double const theta) const { return EvaluateIsotropicKernel(theta); }
