#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/UniformBox.h>
#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/NumericalSolvers/GraphLaplacian/GraphLaplacian.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

TEST(GraphLaplacianExamples, ComputeHeatEigenvalues) {
  // the dimension of the problem
  const std::size_t dim = 2;

  // create a standard Gaussian random variable
  std::vector<std::pair<double, double> > bounds(dim, std::pair<double, double>(1.0, 0.0));
  auto rv = std::make_shared<UniformBox>(bounds)->AsVariable();

  // the options for the graph Laplacian
  YAML::Node options;
  options["NumSamples"] = 10000;
  options["Bandwidth"] = 0.25;
  options["EigensolverTol"] = 1.0e-4;

  // set the kernel options
  YAML::Node kernelOptions;
  kernelOptions["Kernel"] = "HatKernel";
  options["KernelOptions"] = kernelOptions;

  // create the graph laplacian
  auto laplacian = std::make_shared<GraphLaplacian>(rv, options);

  // write the samples to file
  laplacian->WriteToFile("output.h5");

  // construct the heat matrix
  laplacian->ConstructHeatMatrix();

  // compute the largest eigenvalues
  const std::size_t neig = 10;
  const Eigen::VectorXd eigenvalues = laplacian->HeatMatrixEigenvalues(neig);

  std::cout << eigenvalues.transpose() << std::endl;

  // the largest eigenvalue is 1
  EXPECT_NEAR(eigenvalues(0), 1.0, 10.0*laplacian->EigensolverTolerance());
}
