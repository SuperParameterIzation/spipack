#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/UniformBox.h>
#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/NumericalSolvers/GraphLaplacian/GraphLaplacian.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

TEST(GraphLaplacianExamples, ComputeHeatEigenvalues) {
  // the dimension of the problem
  const std::size_t dim = 6;

  // create a standard Gaussian random variable
  std::vector<std::pair<double, double> > bounds(dim, std::pair<double, double>(1.0, 0.0));
  auto rv = std::make_shared<UniformBox>(bounds)->AsVariable();

  // the options for the graph Laplacian
  YAML::Node options;
  options["NumSamples"] = 1000;
  options["Bandwidth"] = 0.25;
  //options["EigenSolverMaxIt"] = 10000;
  options["EigenSolverTol"] = 1.0e-4;

  // set the kernel options
  YAML::Node kernelOptions;
  kernelOptions["Kernel"] = "HatKernel";
  options["KernelOptions"] = kernelOptions;

  // create the graph laplacian
  auto laplacian = std::make_shared<GraphLaplacian>(rv, options);
}
