#include "spipack/Tools/Kernels/IsotropicKernel.hpp"

#include <iostream>

using namespace spi::Tools;

IsotropicKernel::IsotropicKernel(YAML::Node const& options) : Kernel(options) {}

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
