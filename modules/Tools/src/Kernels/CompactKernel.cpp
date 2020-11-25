#include "spipack/Tools/Kernels/CompactKernel.hpp"

#include <iostream>

#include <MUQ/Approximation/Quadrature/GaussQuadrature.h>
#include <MUQ/Approximation/Quadrature/FullTensorQuadrature.h>
#include <MUQ/Approximation/Polynomials/Legendre.h>


using namespace muq::Approximation;
using namespace spi::Tools;

CompactKernel::CompactKernel(YAML::Node const& options) : IsotropicKernel(options) {}

std::shared_ptr<CompactKernel> CompactKernel::Construct(YAML::Node const& options) {
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

std::shared_ptr<CompactKernel::ConstructKernelMap> CompactKernel::KernelMap() {
  // define a static map from type to constructor
  static std::shared_ptr<ConstructKernelMap> map;

  if( !map ) { // if the map has not yet been created ...
    // ... create the map
    map = std::make_shared<ConstructKernelMap>();
  }

  return map;
}

double CompactKernel::EvaluateIsotropicKernel(double const theta) const { return EvaluateCompactKernel(theta); }

double CompactKernel::EvaluateCompactKernel(double const theta) const {
  static const double buffer = 1.0+1.0e-10;
  return (theta<buffer? EvaluateCompactKernelImpl(theta) : 0.0);
}

double CompactKernel::NumericallyIntegrate(std::size_t const dim, std::size_t const n) const {
  // the 1d quadrature rule
  auto poly = std::make_shared<Legendre>();
  auto gq = std::make_shared<GaussQuadrature>(poly);

  // comptue the full tensor quadrature
  FullTensorQuadrature quad(dim, gq, order);

  // get the points and weights
  const Eigen::MatrixXd quadPts = quad.Points();
  const Eigen::VectorXd quadWts = quad.Weights();
  assert(quadPts.cols()==quadWts.size());
  assert(quadPts.rows()==dim);

  // compute the integral
  double integral = 0.0;
  for( std::size_t i=0; i<quadWts.size(); ++i ) {
    const Eigen::VectorXd pt = (quadPts.col(i)+Eigen::VectorXd::Ones(dim))/2.0;

    double pw = 1.0;
    for( std::size_t j=0; j<n; ++j ) { pw *= pt(0); }

    integral += quadWts(i)*pw*EvaluateIsotropicKernel(pt.dot(pt));
  }

  return integral/2.0;
}
