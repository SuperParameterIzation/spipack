#include "spipack/Tools/Kernels/Kernel.hpp"

using namespace spi::Tools;

Kernel::Kernel(YAML::Node const& options) {}

double Kernel::operator()(Eigen::Ref<const Eigen::VectorXd> const& x1, Eigen::Ref<const Eigen::VectorXd> const& x2) const {
  assert(x1.size()==x2.size());
  return Evaluate(x1, x2);
}

std::shared_ptr<Kernel> Kernel::Construct(YAML::Node const& options) {
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

std::shared_ptr<Kernel::ConstructKernelMap> Kernel::KernelMap() {
  // define a static map from type to constructor
  static std::shared_ptr<ConstructKernelMap> map;

  if( !map ) { // if the map has not yet been created ...
    // ... create the map
    map = std::make_shared<ConstructKernelMap>();
  }

  return map;
}
