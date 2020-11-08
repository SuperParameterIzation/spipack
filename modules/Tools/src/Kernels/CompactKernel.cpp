#include "spipack/Tools/Kernels/CompactKernel.hpp"

#include <iostream>

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
