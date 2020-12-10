#include "spipack/Tools/Kernels/IsotropicKernel.hpp"

#include <iostream>

#include <MUQ/Approximation/Quadrature/GaussQuadrature.h>
#include <MUQ/Approximation/Quadrature/FullTensorQuadrature.h>
#include <MUQ/Approximation/Polynomials/PhysicistHermite.h>

#include "spipack/Tools/Kernels/ExponentialKernel.hpp"

using namespace muq::Approximation;
using namespace spi::Tools;

IsotropicKernel::IsotropicKernel(YAML::Node const& options) :
Kernel(options),
order(options["QuadratureOrder"].as<std::size_t>(5)),
delta(options["FiniteDifferenceParameter"].as<double>(1.0e-6))
{}

std::shared_ptr<IsotropicKernel> IsotropicKernel::Construct(YAML::Node const& options) {
  // get the name of the kernel
  const std::string& kernelName = options["Kernel"].as<std::string>();

  // construct it from the map
  auto kernelMap = KernelMap();
  auto iter = kernelMap->find(kernelName);
  if( iter==kernelMap->end() ){
    std::cerr << "ERROR: Could not find spi::Tools::Kernel \"" << kernelName << "\", available types are:\n";

    for(auto it=kernelMap->begin(); it!=kernelMap->end(); ++it)
      std::cerr << "  " << it->first << std::endl;
    std::cerr << std::endl;

    assert(iter != kernelMap->end());
  }

  return iter->second(options);
}


std::shared_ptr<IsotropicKernel::ConstructKernelMap> IsotropicKernel::KernelMap() {
  // define a static map from type to constructor
  static std::shared_ptr<ConstructKernelMap> map;

  if( !map ) { // if the map has not yet been created ...
    // ... create the map
    map = std::make_shared<ConstructKernelMap>();
  }

  return map;
}

double IsotropicKernel::Evaluate(Eigen::Ref<const Eigen::VectorXd> const& x1, Eigen::Ref<const Eigen::VectorXd> const& x2) const {
  assert(x1.size()==x2.size());
  const Eigen::VectorXd diff = x1-x2;
  return EvaluateIsotropicKernel(diff.dot(diff));
}

double IsotropicKernel::operator()(double const theta) const { return EvaluateIsotropicKernel(theta); }

double IsotropicKernel::IsotropicKernelDerivative(double const theta) const { return IsotropicKernelDerivativeFD(theta); }

double IsotropicKernel::IsotropicKernelSecondDerivative(double const theta) const { return IsotropicKernelSecondDerivativeFD(theta); }

double IsotropicKernel::IsotropicKernelDerivativeFD(double const theta) const { return (EvaluateIsotropicKernel(theta+delta)-EvaluateIsotropicKernel(theta))/delta; }

double IsotropicKernel::IsotropicKernelSecondDerivativeFD(double const theta) const {
  return (IsotropicKernelDerivative(theta+delta)-IsotropicKernelDerivative(theta))/delta;
}

double IsotropicKernel::Integrate(std::size_t const dim, std::size_t const n) const { return NumericallyIntegrate(dim, n); }

double IsotropicKernel::NumericallyIntegrate(std::size_t const dim, std::size_t const n) const {
  // the 1d quadrature rule
  auto poly = std::make_shared<PhysicistHermite>();
  auto gq = std::make_shared<GaussQuadrature>(poly);

  // comptue the full tensor quadrature
  FullTensorQuadrature quad(dim, gq, order);

  // get the points and weights
  const Eigen::MatrixXd quadPts = quad.Points();
  const Eigen::VectorXd quadWts = quad.Weights();
  assert(quadPts.cols()==quadWts.size());

  // need to normalize the kernel since we are using Guass-Hermite
  YAML::Node options;
  const ExponentialKernel expkern(options);

  // compute the integral
  double integral = 0.0;
  for( std::size_t i=0; i<quadWts.size(); ++i ) {
    double pw = 1.0;
    for( std::size_t j=0; j<n; ++j ) { pw *= quadPts(0, i); }

    const double theta = quadPts.col(i).dot(quadPts.col(i));
    integral += quadWts(i)*pw*EvaluateIsotropicKernel(theta)/expkern.EvaluateIsotropicKernel(theta);
  }

  return integral;
}
